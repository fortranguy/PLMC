module mod_physique

use mod_tools
use data_cell
use data_particles
use data_mc
use data_potentiel
use data_constants

implicit none

    contains
    
    ! Condition initiale ------------------------------------------------------
    
    subroutine condIni(unitRapport)
    
        integer, intent(in) :: unitRapport
        
        real(DP) :: compac, densite
        character(len=20) :: init
        integer :: longueur, statut
        
        call get_command_argument(1, init, longueur, statut)
		if (statut /= 0) stop "erreur get_command_argument"
		if (command_argument_count() > 1) stop "Trop d'arguments"
		    
        write(unitRapport, *) "Condition initiale :"
        
        select case (init)
            case ("cube")
                call iniPosCub()
                write(unitRapport, *) "    Cubique simple"
            case ("alea")
                call iniPosAlea()
                write(unitRapport, *) "    Déposition aléatoire"
            case default
                write(*, *) "Préciser la condition initiale : "
                write(*, *) "   'cube' ou 'alea'."
                stop
        end select
        
        densite = real(Ncol1, DP) / product(Lsize)
        write(*, *) "    Densité = ", densite
        write(unitRapport, *) "    Densité = ", densite
        
        compac = 4._DP/3._DP*PI*rayon1**3 * Ncol1 / product(Lsize)
        write(*, *) "    Compacité = ", compac
        write(unitRapport, *) "    Compacité = ", compac
        
    end subroutine condIni
    
    subroutine iniPosCub()
    
        integer :: iDir
        integer :: i, j, k, iPart
        integer, dimension(Dim) :: nParts
        real(DP), dimension(Dim) :: ratio
        real(DP) :: oneThird = 1._DP/3._DP
        
        write(*, *) "Cubique simple"
        
        ! Proportion selon la direction
        
        nParts(1) = int( (Ncol1*Lsize(1)**2/Lsize(2)/Lsize(3))**oneThird )
        nParts(2) = int( (Ncol1*Lsize(2)**2/Lsize(3)/Lsize(1))**oneThird )
        nParts(3) = int( (Ncol1*Lsize(3)**2/Lsize(1)/Lsize(2))**oneThird )
        
        ! Vérification
        
        iDir = 1
        do while (product(nParts)<Ncol1)
        	nParts(iDir) = nParts(iDir) + 1
        	iDir = iDir + 1
        end do
        
        ratio(:) = Lsize(:)/real(nParts(:), DP) ! A vérifier
        do iDir = 1, Dim
            if ( rmin*real(nParts(iDir), DP) > Lsize(iDir) ) then
                write(*, *) "    Problème : trop dense dans la direction ",&
                iDir
                stop
            end if
        end do
        
        ! Remplissage
        
        do k = 1, nParts(3)
            do j = 1, nParts(2)
                do i = 1, nParts(1)
                    if (i*j*k <= Ncol1) then
                        iPart = i + nParts(1)*(j-1) + nParts(1)*nParts(2)*(k-1)
                        X(1, iPart) = ratio(1)*real(i, DP)
                        X(2, iPart) = ratio(2)*real(j, DP)
                        X(3, iPart) = ratio(3)*real(k, DP)
                        ! A vérifier
                    end if
                end do
            end do
        end do
    
    	do iDir = 1, Dim
            x(iDir, :) = x(iDir, :) - 0.5_DP*ratio(iDir) ! nécessaire ? tradition
        end do
    
    end subroutine iniPosCub
    
    ! ---------------------------------
    
    subroutine iniPosAlea()
    
        integer :: iPart, Nparts, nOK
        real(DP), dimension(Dim) :: xTest, DeltaX
    
        write(*, *) "Déposition aléatoire"
    
        call random_number(X(:, 1))
        X(:, 1) = X(:, 1)*(Lsize(:)-2*Rayon1)
        Nparts = 1        
        
        do while (Nparts<Ncol1)
        
            call random_number(xTest)
            xTest(:) = xTest(:)*(Lsize(:)-2._DP*Rayon1)
            
            nOK = 0
            do iPart = 1, Nparts
                DeltaX(:) = X(:, iPart) - xTest(:)
                call pbc_dif(DeltaX)
                if( sqrt(dot_product(DeltaX, DeltaX)) >= rmin) then
                    nOK = nOK + 1
                else
                    exit
                end if
            end do
            
            if (nOK == Nparts) then
                Nparts = Nparts + 1
                X(:, Nparts) = xTest(:)
                write(*, *) "    Particule", Nparts, "OK"
            end if
            
        end do
        
        do iPart = 1, Ncol1
            X(:, iPart) = X(:, iPart) + Rayon1
        end do
    
    end subroutine iniPosAlea
    
    ! Test d'overlap ----------------------------------------------------------
    
    subroutine overlapTest()
    
        integer :: jPart, iPart
        real(DP), dimension(Dim) :: DeltaX
    
        do jPart = 1, Ncol1
            do iPart = 1, Ncol1
                if (iPart /= jPart) then
                    
                    DeltaX(:) = X(:, iPart) - X(:, jPart)
                    call pbc_dif(DeltaX)
                    
                    if (sqrt(dot_product(DeltaX, DeltaX)) < rmin) then
                        write(*, *) "    Overlap !", iPart, jPart
                        write(*, * ) "    r_ij = ", &
                        sqrt(dot_product(DeltaX, DeltaX))
                        stop
                    end if
                    
                end if
            end do
        end do
        
        write(*, *) "    Overlap test : OK !"
    
    end subroutine overlapTest
    
    ! Conditions aux limites périodiques (PBC) --------------------------------
    
    subroutine pbc_dif(DeltaX)
    
        real(DP), dimension(Dim), intent(inout) :: DeltaX
        
        DeltaX(:) = modulo(DeltaX(:), Lsize(:))
        
        where( DeltaX(:) > LsizeMi(:) )
            DeltaX(:) = DeltaX(:) - Lsize(:)
        end where
    
    end subroutine pbc_dif
    
    ! Energie potentielle -----------------------------------------------------
    
    subroutine ePotIni()
    
        integer :: i
        real(DP) :: r
        
        do i = iMin, Ntab11
        
            r = rcut11*real(i, DP)/real(Ntab11, DP)
            ! Yukawa :
            Vtab11(i) = epsilon11*exp(-alpha11*(r-rmin))/r
            
        end do
    
    end subroutine ePotIni
    
    ! ---------------------------------
    
    function ePot(r)
    
        real(DP), intent(in) :: r
        
        integer :: i
        real(DP) :: ri, riPlus
        real(DP) :: ePot
        
        if (r < rcut11) then
        
            i = int(r/rcut11*real(Ntab11, DP))
            ri = rcut11*real(i, DP)/real(Ntab11, DP)
            riPlus = rcut11*real(i+1, DP)/real(Ntab11, DP)
            
            ePot = Vtab11(i) + (r-ri)/(riPlus-ri)*(Vtab11(i+1)-Vtab11(i))
            
        else
        
            ePot = 0._DP
            
        end if
    
    end function ePot    
    
    ! ---------------------------------
    
    subroutine ePotDif(iOld, xNew, overlap, dEn)
    
        integer, intent(in) :: iOld
        real(DP), dimension(Dim), intent(in) :: xNew
        logical, intent(inout) :: overlap
        real(DP), intent(inout) :: dEn
        
        integer :: jPart
        real(DP), dimension(Dim) :: DeltaXnew, DeltaXold
        real(DP) rNew, rOld
        
        overlap = .false.
        dEn = 0._DP
        
        do jPart = 1, Ncol1
            if (jPart /= iOld) then
            
                DeltaXnew(:) = X(:, jPart) - xNew(:)
                call pbc_dif(DeltaXnew)
                DeltaXold(:) = X(:, jPart) - X(:, iOld)
                call pbc_dif(DeltaXold)
                
                rNew = sqrt(dot_product(DeltaXnew, DeltaXnew))
                if (rNew < rmin) then
                    overlap = .true.
                    return
                end if                
                rOld = sqrt(dot_product(DeltaXold, DeltaXold))
                
                dEn = dEn + ePot(rNew) - ePot(rOld)    
            
            end if
        end do
    
    end subroutine ePotDif
    
    ! Déplacement d'une particule ---------------------------------------------
    
    subroutine mcMove(enTot, Nrejects)
    
        real(DP), intent(inout) :: enTot
        integer, intent(inout) :: Nrejects
        
        logical :: overlap
        integer :: iOld
        real(DP) :: rand
        real(DP), dimension(Dim) :: xNew
        real(DP) :: dEn
        
        call random_number(rand)
        iOld = int(rand*Ncol1) + 1
        call random_number(xNew)
        xNew(:) = X(:, iOld) + (xNew(:)-0.5_DP)*dx(:)
        xNew(:) = modulo(xNew(:), Lsize(:))

        call ePotDif(iOld, xNew, overlap, dEn)
        
        if (.not. overlap) then
            call random_number(rand)
            if ( rand < exp(-dEn/Tstar) ) then
                X(:, iOld) = xNew(:)
                enTot = enTot + dEn
            else
                Nrejects = Nrejects + 1
            end if
        else
            Nrejects = Nrejects + 1
        end if
    
    end subroutine mcMove
    
    ! Méthode de Widom --------------------------------------------------------
    
    subroutine ePotTest(xTest, overlap, enTest)

        real(DP), dimension(Dim), intent(in) :: xTest
        real(DP), intent(inout) :: enTest
        logical, intent(inout) :: overlap
        
        integer :: iPart
        real(DP), dimension(Dim) :: DeltaXtest
        real(DP) :: rTest
        
        overlap = .false.
        enTest = 0._DP
        
            do iPart = 1, Ncol1
                
                DeltaXtest(:) = X(:, iPart) - xTest(:)
                call pbc_dif(DeltaXtest)
                        
                rTest = sqrt(dot_product(DeltaXtest, DeltaXtest))
                if (rTest < rmin) then
                    overlap = .true.
                    return
                end if
                
                enTest = enTest + ePot(rTest)   
                
            end do
        
    end subroutine ePotTest
        
    ! ---------------------------------
    
    subroutine widom(nWidom, Lratio, activExInv)
        
        integer, intent(in) :: nWidom
        real(DP), intent(in) :: Lratio
        real(DP), intent(inOut) :: activExInv 
        
        integer :: iWid
        real(DP) :: widTestSum
        real(DP), dimension(Dim) :: xTest
        logical :: overlap        
        real(DP) :: enTest
        
        widTestSum = 0._DP
        
        do iWid = 1, nWidom           
            
            call random_number(xTest)
            xTest(:) = Lratio*Lsize(:) * xTest(:)            
            call ePotTest(xTest, overlap, enTest)
            
            if (.not. overlap) then
                widTestSum = widTestSum + exp(-enTest/Tstar)
            end if
            
        end do
        
        activExInv = widTestSum/real(nWidom, DP)
        
    end subroutine widom
    
    ! Test de consistance -----------------------------------------------------
    
    function enTotCalc()
    
        integer :: iPart, jPart
        real(DP) :: r_ij
        real(DP), dimension(Dim) :: DeltaX
        real(DP) :: enTotCalc
    
        enTotCalc = 0._DP
        
        do jPart = 1, Ncol1
            do iPart = 1, Ncol1
                if (iPart /= jPart) then
                
                    DeltaX(:) = X(:, jPart) - X(:, iPart)
                    call pbc_dif(DeltaX)
                    r_ij = sqrt(dot_product(DeltaX, DeltaX))
                    
                    enTotCalc = enTotCalc + ePot(r_ij)
                    
                end if
            end do
        end do
        
        enTotCalc = 0.5_DP*enTotCalc
    
    end function enTotCalc
    
    ! ---------------------------------
    
    subroutine consisTest(enTot, unitRapport)
    
        real(DP), intent(in) :: enTot
        integer, intent(in) :: unitRapport
        
        write(unitRapport, *) "Test de consistence :"
        write(unitRapport, *) "    enTot_mc_c = ", enTot
        write(unitRapport, *) "    enTot_calc = ", enTotCalc()
    
    end subroutine consisTest

end module mod_physique
