module mod_tools

use data_constants
use data_particles
use data_mc
use data_potentiel
use data_neighbours
use mod_pbc
use class_component

implicit none

    contains

    ! Générateurs de nombres aléatoires : graine ------------------------------
    
    subroutine init_random_seed(unitRapport)
    
        integer, intent(in) :: unitRapport
    
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed(:) = clock + 37 * [ (i - 1, i = 1, n) ]
        call random_seed(put = seed)
        
        write(unitRapport, *) "Graine :"
        write(unitRapport ,*) "    n = ", n
        write(unitRapport ,*) "    seed(:) = ", seed(:)

        deallocate(seed)
        
    end subroutine init_random_seed
    
    ! Condition initiale ------------------------------------------------------
    
    subroutine condIni(unitRapport, sph_X)
    
        integer, intent(in) :: unitRapport
        real(DP), dimension(:, :), intent(inout) :: sph_X
        
        real(DP) :: compac, densite
        character(len=20) :: init
        integer :: longueur, statut
        
        call get_command_argument(1, init, longueur, statut)
        if (statut /= 0) stop "erreur get_command_argument"
        if (command_argument_count() > 1) stop "Trop d'arguments"
            
        write(unitRapport, *) "Condition initiale :"
        
        select case (init)
            case ("cube")
                call iniPosCub(sph_X)
                write(unitRapport, *) "    Cubique simple"
            case ("alea")
                call iniPosAlea(sph_X)
                write(unitRapport, *) "    Déposition aléatoire"
            case default
                write(*, *) "Préciser la condition initiale : "
                write(*, *) "   'cube' ou 'alea'."
                stop
        end select
        
        densite = real(sph_Ncol, DP) / product(Lsize)
        write(*, *) "    Densité = ", densite
        write(unitRapport, *) "    Densité = ", densite
        
        compac = 4._DP/3._DP*PI*sph_radius**3 * densite
        write(*, *) "    Compacité = ", compac
        write(unitRapport, *) "    Compacité = ", compac
        
    end subroutine condIni
    
    subroutine iniPosCub(sph_X)
    
        real(DP), dimension(:, :), intent(inout) :: sph_X
    
        integer :: iDir
        integer :: i, j, k, iCol
        integer, dimension(Dim) :: nCols
        real(DP), dimension(Dim) :: ratio
        real(DP) :: oneThird = 1._DP/3._DP
        
        write(*, *) "Cubique simple"
        
        ! Proportion selon la direction
        
        nCols(1) = int( (sph_Ncol*Lsize(1)**2/Lsize(2)/Lsize(3))**oneThird )
        nCols(2) = int( (sph_Ncol*Lsize(2)**2/Lsize(3)/Lsize(1))**oneThird )
        nCols(3) = int( (sph_Ncol*Lsize(3)**2/Lsize(1)/Lsize(2))**oneThird )
        
        ! Vérification
        
        iDir = 1
        do while (product(nCols)<sph_Ncol)
            nCols(iDir) = nCols(iDir) + 1
            iDir = iDir + 1
        end do
        
        ratio(:) = Lsize(:)/real(nCols(:), DP) ! A vérifier
        do iDir = 1, Dim
            if ( sph_rmin*real(nCols(iDir), DP) > Lsize(iDir) ) then
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
                    if (iCol <= sph_Ncol) then
                        sph_X(1, iCol) = ratio(1)*real(i, DP)
                        sph_X(2, iCol) = ratio(2)*real(j, DP)
                        sph_X(3, iCol) = ratio(3)*real(k, DP)
                    end if
                end do
            end do
        end do
    
        do iDir = 1, Dim
            sph_X(iDir, :) = sph_X(iDir, :) - 0.5_DP*ratio(iDir) ! nécessaire ? tradition
        end do
    
    end subroutine iniPosCub
    
    ! ---------------------------------
    
    subroutine iniPosAlea(sph_X)
    
        real(DP), dimension(:, :), intent(inout) :: sph_X
    
        integer :: iCol, Ncols, nOK
        real(DP), dimension(Dim) :: xTest
        real(DP) :: rTest
    
        write(*, *) "Déposition aléatoire"
    
        call random_number(sph_X(:, 1))
        sph_X(:, 1) = sph_X(:, 1)*(Lsize(:)-2*sph_radius)
        Ncols = 1        
        
        do while (Ncols<sph_Ncol)
        
            call random_number(xTest)
            xTest(:) = xTest(:)*(Lsize(:)-2._DP*sph_radius)
            
            nOK = 0
            do iCol = 1, Ncols
                rTest = dist(sph_X(:, iCol), xTest(:))
                if (rTest >= sph_rmin) then
                    nOK = nOK + 1
                else
                    exit
                end if
            end do
            
            if (nOK == Ncols) then
                Ncols = Ncols + 1
                sph_X(:, Ncols) = xTest(:)
                write(*, *) "    Particule", Ncols, "OK"
            end if
            
        end do
        
        do iCol = 1, sph_Ncol
            sph_X(:, iCol) = sph_X(:, iCol) + sph_radius
        end do
    
    end subroutine iniPosAlea
    
    ! Résultats ---------------------------------------------------------------
        
    subroutine mcResults(enTotSum, activExInvSum, tauxRejectsSum, duree,&
    	unitRapport)

        real(DP), intent(in) :: enTotSum, activExInvSum     
        real(DP), intent(in) :: tauxRejectsSum
        real(DP), intent(in) :: duree
        integer, intent(in) :: unitRapport
        
        real(DP) :: realNstep = real(Nstep, DP)
        real(DP) :: potChiId, potChiEx
    
        write(unitRapport, *) "Résultats :"
        write(unitRapport, *) "    Energie moyenne = ", &
            enTotSum/realNstep
        write(unitRapport, *) "    Energie moyenne par particule = ", &
            enTotSum/realNstep/real(sph_Ncol, DP)
        potChiId = -Tstar*log( product(Lsize)/real(sph_Ncol+1,DP) )
        write(unitRapport, *) "    Potentiel chimique idéal = ", potChiId
        potChiEx = -Tstar*log( activExInvSum/realNstep )
        write(unitRapport, *) "    Potentiel chimique (excès) moyen = ", &
            potChiEx           
        write(unitRapport, *) "    potChi.moy = ", potChiId + potChiEx
        write(unitRapport, *) "    Taux rejets = ", &
            tauxRejectsSum/real(Nstep+Ntherm, DP)
        write(unitRapport, *) "    Durée =", duree/60._DP, "min"        
    
    end subroutine mcResults
    
end module mod_tools
    
