module mod_physique

use data_constants
use data_cell
use data_particles
use data_mc
use data_potentiel
use data_neighbours
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
                stop
        end select
        
        densite = real(Ncol, DP) / product(Lsize)
        write(*, *) "    Densité = ", densite
        write(unitRapport, *) "    Densité = ", densite
        
        compac = 4._DP/3._DP*PI*radius**3 * densite
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
        
        nCols(1) = int( (Ncol*Lsize(1)**2/Lsize(2)/Lsize(3))**oneThird )
        nCols(2) = int( (Ncol*Lsize(2)**2/Lsize(3)/Lsize(1))**oneThird )
        nCols(3) = int( (Ncol*Lsize(3)**2/Lsize(1)/Lsize(2))**oneThird )
        
        ! Vérification
        
        iDir = 1
        do while (product(nCols)<Ncol)
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
                    iCol = i + nCols(1)*(j-1) + nCols(1)*nCols(2)*(k-1)
                    if (iCol <= Ncol) then
                        X(1, iCol) = ratio(1)*real(i, DP)
                        X(2, iCol) = ratio(2)*real(j, DP)
                        X(3, iCol) = ratio(3)*real(k, DP)
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
        real(DP), dimension(Dim) :: xTest
        real(DP) :: rTest
    
        write(*, *) "Déposition aléatoire"
    
        call random_number(X(:, 1))
        X(:, 1) = X(:, 1)*(Lsize(:)-2*radius)
        Ncols = 1        
        
        do while (Ncols<Ncol)
        
            call random_number(xTest)
            xTest(:) = xTest(:)*(Lsize(:)-2._DP*radius)
            
            nOK = 0
            do iCol = 1, Ncols
                rTest = dist(X(:, iCol), xTest(:))
                if (rTest >= rmin) then
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
        
        do iCol = 1, Ncol
            X(:, iCol) = X(:, iCol) + radius
        end do
    
    end subroutine iniPosAlea
    
    ! Test d'overlap ----------------------------------------------------------
    
    subroutine overlapTest()
    
        integer :: jCol, iCol
        real(DP) :: r_ij
    
        do jCol = 1, Ncol
            do iCol = 1, Ncol
                if (iCol /= jCol) then
                    
                    r_ij = dist(X(:, iCol), X(:, jCol))
                    if (r_ij < rmin) then
                        write(*, *) "    Overlap !", iCol, jCol
                        write(*, * ) "    r_ij = ", r_ij
                        stop
                    end if
                    
                end if
            end do
        end do
        
        write(*, *) "    Overlap test : OK !"
    
    end subroutine overlapTest
    
    ! Distance entre 2 particules (CLP) ---------------------------------------
    
    function dist(X1, X2)
    
        real(DP), dimension(Dim), intent(in) :: X1, X2
        real(DP), dimension(Dim) :: DeltaX
        real(DP) :: dist
        
        DeltaX(:) = X2(:) - X1(:)
        DeltaX(:) = modulo(DeltaX(:), Lsize(:))
        
        where( DeltaX(:) > LsizeMi(:) )
            DeltaX(:) = DeltaX(:) - Lsize(:)
        end where
        
        dist = sqrt(dot_product(DeltaX, DeltaX))
    
    end function dist

end module mod_physique
