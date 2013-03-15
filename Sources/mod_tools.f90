module mod_tools

use data_constants
use data_particles
use data_mc
use mod_physics

implicit none

contains

    ! Random number generator : seed ------------------------------
    
    subroutine init_random_seed(unitRapport)
    
        integer, intent(in) :: unitRapport
    
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed(:) = clock + 37 * [ (i - 1, i = 1, n) ]
        call random_seed(put = seed)
        
        write(unitRapport, *) "RNG :"
        write(unitRapport ,*) "    n = ", n
        write(unitRapport ,*) "    seed(:) = ", seed(:)

        deallocate(seed)
        
    end subroutine init_random_seed
    
    ! Initial condition ------------------------------------------------------
    
    subroutine condIni(unitRapport, sph_X)
    
        integer, intent(in) :: unitRapport
        real(DP), dimension(:, :), intent(inout) :: sph_X
        
        real(DP) :: compac, densite
        character(len=20) :: init
        integer :: longueur, statut
        
        call get_command_argument(1, init, longueur, statut)
        if (statut /= 0) stop "error get_command_argument"
        if (command_argument_count() > 1) stop "Too many arguments"
            
        write(unitRapport, *) "Initial condition :"
        
        select case (init)
            case ("cube")
                call iniPosCub(sph_X)
                write(unitRapport, *) "    Primitive cubic"
            case ("rand")
                call iniPosAlea(sph_X)
                write(unitRapport, *) "    Random deposition"
            case default
                write(*, *) "Enter the initial condition : "
                write(*, *) "   'cube' or 'rand'."
                stop
        end select
        
        densite = real(sph_Ncol, DP) / product(Lsize)
        write(*, *) "    Density = ", densite
        write(unitRapport, *) "    Density = ", densite
        
        compac = 4._DP/3._DP*PI*sph_radius**3 * densite
        write(*, *) "    Compacity = ", compac
        write(unitRapport, *) "    Compacity = ", compac
        
    end subroutine condIni
    
    subroutine iniPosCub(sph_X)
    
        real(DP), dimension(:, :), intent(inout) :: sph_X
    
        integer :: iDir
        integer :: i, j, k, iCol
        integer, dimension(Dim) :: nCols
        real(DP), dimension(Dim) :: ratio
        real(DP) :: oneThird = 1._DP/3._DP
        
        write(*, *) "Primitive cubic"
        
        ! Proportion according to the direction
        
        nCols(1) = int( (sph_Ncol*Lsize(1)**2/Lsize(2)/Lsize(3))**oneThird )
        nCols(2) = int( (sph_Ncol*Lsize(2)**2/Lsize(3)/Lsize(1))**oneThird )
        nCols(3) = int( (sph_Ncol*Lsize(3)**2/Lsize(1)/Lsize(2))**oneThird )
        
        ! Check
        
        iDir = 1
        do while (product(nCols)<sph_Ncol)
            nCols(iDir) = nCols(iDir) + 1
            iDir = iDir + 1
        end do
        
        ratio(:) = Lsize(:)/real(nCols(:), DP) ! A vérifier
        do iDir = 1, Dim
            if ( sph_rmin*real(nCols(iDir), DP) > Lsize(iDir) ) then
                write(*, *) "    Error : too dense in the direction", iDir
                stop
            end if
        end do
        
        ! Filling
        
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
            sph_X(iDir, :) = sph_X(iDir, :) - 0.5_DP*ratio(iDir) ! just inside
        end do
    
    end subroutine iniPosCub
    
    ! ---------------------------------
    
    subroutine iniPosAlea(sph_X)
    
        real(DP), dimension(:, :), intent(inout) :: sph_X
    
        integer :: iCol, Ncols, nOK
        real(DP), dimension(Dim) :: xTest
        real(DP) :: rTest
    
        write(*, *) "Random deposition"
    
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
    
    ! Results ---------------------------------------------------------------
        
    subroutine mcResults(enTotSum, activExInvSum, tauxRejectsSum, duration, &
    	unitRapport)

        real(DP), intent(in) :: enTotSum, activExInvSum     
        real(DP), intent(in) :: tauxRejectsSum
        real(DP), intent(in) :: duration
        integer, intent(in) :: unitRapport
        
        real(DP) :: realNstep = real(Nstep, DP)
        real(DP) :: potChiId, potChiEx
    
        write(unitRapport, *) "Results :"
        write(unitRapport, *) "    average energy = ", &
            enTotSum/realNstep
        write(unitRapport, *) "    average energy per particule = ", &
            enTotSum/realNstep/real(sph_Ncol, DP)
        potChiId = -Tstar*log( product(Lsize)/real(sph_Ncol+1,DP) )
        write(unitRapport, *) "    ideal chemical potential = ", potChiId
        potChiEx = -Tstar*log( activExInvSum/realNstep )
        write(unitRapport, *) "    average excess chemical potential = ", &
            potChiEx           
        write(unitRapport, *) "    potChi.avg = ", potChiId + potChiEx
        write(unitRapport, *) "    Rejection rate = ", &
            tauxRejectsSum/real(Nstep+Ntherm, DP)
        write(unitRapport, *) "    duration =", duration/60._DP, "min"        
    
    end subroutine mcResults
    
end module mod_tools
