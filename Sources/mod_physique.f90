module mod_physique

use data_constants
use data_cell
use data_particles
use data_mc
use data_potentiel
use mod_neighbours

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
        end select
        
        densite = real(Ncol1, DP) / product(Lsize)
        write(*, *) "    Densité = ", densite
        write(unitRapport, *) "    Densité = ", densite
        
        compac = 4._DP/3._DP*PI*rayon1**3 * densite
        write(*, *) "    Compacité = ", compac
        write(unitRapport, *) "    Compacité = ", compac
        
    end subroutine condIni
    
    subroutine iniPosCub()
    
        integer :: iDir
        integer :: i, j, k, iCol
        integer, dimension(Dim) :: nCols
        real(DP), dimension(Dim) :: ratio
        real(DP) :: oneThird = 1._DP/3._DP
        
        write(*, *) "Cubique simple"
        
        ! Proportion selon la direction
        
        nCols(1) = int( (Ncol1*Lsize(1)**2/Lsize(2)/Lsize(3))**oneThird )
        nCols(2) = int( (Ncol1*Lsize(2)**2/Lsize(3)/Lsize(1))**oneThird )
        nCols(3) = int( (Ncol1*Lsize(3)**2/Lsize(1)/Lsize(2))**oneThird )
        
        ! Vérification
        
        iDir = 1
        do while (product(nCols)<Ncol1)
            nCols(iDir) = nCols(iDir) + 1
            iDir = iDir + 1
        end do
        
        ratio(:) = Lsize(:)/real(nCols(:), DP) ! A vérifier
        do iDir = 1, Dim
            if ( rmin*real(nCols(iDir), DP) > Lsize(iDir) ) then
                write(*, *) "    Problème : trop dense dans la direction ",&
                iDir
                stop
            end if
        end do
        
        ! Remplissage
        
        do k = 1, nCols(3)
            do j = 1, nCols(2)
                do i = 1, nCols(1)
                    if (i*j*k <= Ncol1) then
                        iCol = i + nCols(1)*(j-1) + nCols(1)*nCols(2)*(k-1)
                        X(1, iCol) = ratio(1)*real(i, DP)
                        X(2, iCol) = ratio(2)*real(j, DP)
                        X(3, iCol) = ratio(3)*real(k, DP)
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
    
        integer :: iCol, Ncols, nOK
        real(DP), dimension(Dim) :: xTest, DeltaX
    
        write(*, *) "Déposition aléatoire"
    
        call random_number(X(:, 1))
        X(:, 1) = X(:, 1)*(Lsize(:)-2*Rayon1)
        Ncols = 1        
        
        do while (Ncols<Ncol1)
        
            call random_number(xTest)
            xTest(:) = xTest(:)*(Lsize(:)-2._DP*Rayon1)
            
            nOK = 0
            do iCol = 1, Ncols
                DeltaX(:) = X(:, iCol) - xTest(:)
                call pbc_dif(DeltaX)
                if( sqrt(dot_product(DeltaX, DeltaX)) >= rmin) then
                    nOK = nOK + 1
                else
                    exit
                end if
            end do
            
            if (nOK == Ncols) then
                Ncols = Ncols + 1
                X(:, Ncols) = xTest(:)
                write(*, *) "    Particule", Ncols, "OK"
            end if
            
        end do
        
        do iCol = 1, Ncol1
            X(:, iCol) = X(:, iCol) + Rayon1
        end do
    
    end subroutine iniPosAlea
    
    ! Test d'overlap ----------------------------------------------------------
    
    subroutine overlapTest()
    
        integer :: jCol, iCol
        real(DP), dimension(Dim) :: DeltaX
    
        do jCol = 1, Ncol1
            do iCol = 1, Ncol1
                if (iCol /= jCol) then
                    
                    DeltaX(:) = X(:, iCol) - X(:, jCol)
                    call pbc_dif(DeltaX)
                    
                    if (sqrt(dot_product(DeltaX, DeltaX)) < rmin) then
                        write(*, *) "    Overlap !", iCol, jCol
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
   
    function ePot(r)
   
        real(DP), intent(in) :: r
       
        integer :: i
        real(DP) :: r_i, ePot
       
        if (r < rcut11) then
       
            i = int(r/rcut11*real(Ntab11, DP))
            r_i = real(i, DP)*pas11
            ePot = Vtab11(i) + (r-r_i)/pas11 * (Vtab11(i+1)-Vtab11(i))
           
        else
       
            ePot = 0._DP
           
        end if
   
    end function ePot
    
    subroutine ePotNeigh(iCol, pos, iCell, overlap, energ)
        
        integer, intent(in) :: iCol, iCell
        real(DP), dimension(Dim), intent(in) :: pos
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNeigh,  iCell_neigh
        real(DP), dimension(Dim) :: DeltaX
        real(DP) :: r
    
        type(Particle), pointer :: courant => null(), suivant => null()
        
        overlap = .false.
        energ = 0._DP
    
        do iNeigh = 1, cell_neighs_nb
        
            iCell_neigh = cell_neighs(iNeigh, iCell)
            courant => cellsBegin(iCell_neigh)%particle%next            
            if (.not. associated(courant%next)) cycle
            
            do
            
                suivant => courant%next
            
                if (courant%iCol /= iCol) then
                
                    DeltaX(:) = X(:, courant%iCol) - pos(:)
                    call pbc_dif(DeltaX)
                    r = sqrt(dot_product(DeltaX, DeltaX))
                    if (r < rmin) then
                        overlap = .true.
                        return
                    end if
                    energ = energ + ePot(r)
       
                end if
                
                if (.not. associated(suivant%next)) exit
                
                courant => suivant
            
            end do            
            
        end do
    
    end subroutine ePotNeigh
    
    ! Déplacement d'une particule ---------------------------------------------
    
    subroutine mcMove(enTot, Nrejects)
    
        real(DP), intent(inout) :: enTot
        integer, intent(inout) :: Nrejects
        
        logical :: overlap
        integer :: iOld
        real(DP) :: rand
        real(DP), dimension(Dim) :: xNew
        integer :: iCellBefore, iCellAfter
        real(DP) :: eNew, eOld, dEn
        
        call random_number(rand)
        iOld = int(rand*Ncol1) + 1
        
        call random_number(xNew)
        xNew(:) = X(:, iOld) + (xNew(:)-0.5_DP)*dx(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        iCellAfter = position_to_cell(xNew)
        call ePotNeigh(iOld, xNew, iCellAfter, overlap, eNew)
        
        if (.not. overlap) then
        
            iCellBefore = col_to_cell(iOld)
            call ePotNeigh(iOld, x(:, iOld), iCellBefore, overlap, eOld)
        	
            dEn = eNew - eOld
        
            call random_number(rand)
            if ( rand < exp(-dEn/Tstar) ) then
                X(:, iOld) = xNew(:)
                enTot = enTot + dEn
                
                if ( iCellBefore /= iCellAfter ) then                
                    call remove_cell_col(iOld, iCellBefore)
                    call add_cell_col(iOld, iCellAfter)
                end if
                
            else
                Nrejects = Nrejects + 1
            end if
            
        else
        
            Nrejects = Nrejects + 1
            
        end if
    
    end subroutine mcMove
    
    ! Méthode de Widom --------------------------------------------------------

    subroutine widom(nWidom, Lratio, activExInv)
        
        integer, intent(in) :: nWidom
        real(DP), intent(in) :: Lratio
        real(DP), intent(inOut) :: activExInv 
        
        integer :: iWid
        real(DP) :: widTestSum
        real(DP), dimension(Dim) :: xTest
        integer :: iCellTest
        logical :: overlap        
        real(DP) :: enTest
        
        widTestSum = 0._DP
        
        do iWid = 1, nWidom           
            
            call random_number(xTest)
            xTest(:) = Lratio*Lsize(:) * xTest(:)    
            iCellTest = position_to_cell(xTest)
            call ePotNeigh(0, xTest, iCellTest, overlap, enTest) 
            
            if (.not. overlap) then
                widTestSum = widTestSum + exp(-enTest/Tstar)
            end if
            
        end do
        
        activExInv = widTestSum/real(nWidom, DP)
        
    end subroutine widom
    
    ! Test de consistance -----------------------------------------------------
    
    function enTotCalc()
    
        integer :: iCol, jCol
        real(DP) :: r_ij
        real(DP), dimension(Dim) :: DeltaX
        real(DP) :: enTotCalc
    
        enTotCalc = 0._DP
        
        do jCol = 1, Ncol1
            do iCol = 1, Ncol1
                if (iCol /= jCol) then
                
                    DeltaX(:) = X(:, jCol) - X(:, iCol)
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
