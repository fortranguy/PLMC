!> \brief Subroutines as tools for the main program

module mod_tools

use data_constants
use data_particles
use data_mc
use mod_physics

implicit none

contains

    !> Random number generator : seed
    
    subroutine init_random_seed(unitReport)
    
        integer, intent(in) :: unitReport
    
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed(:) = clock + 37 * [ (i - 1, i = 1, n) ]
        call random_seed(put = seed)
        
        write(unitReport, *) "RNG :"
        write(unitReport ,*) "    n = ", n
        write(unitReport ,*) "    seed(:) = ", seed(:)

        deallocate(seed)
        
    end subroutine init_random_seed
    
    !> Initial condition
    
    subroutine initialCondition(type1_X, type1_rMin, type2_X, type2_rMin, 
        mix_rMin, unitReport)
    
        real(DP), dimension(:, :), intent(inout) :: type1_X, type2_X
        real(DP), intent(in) :: type1_rMin, type2_rMin, mix_rMin
        integer, intent(in) :: unitReport
        

        character(len=20) :: init
        integer :: longueur, statut
        
        call get_command_argument(1, init, longueur, statut)
        if (statut /= 0) stop "error get_command_argument"
        if (command_argument_count() > 1) stop "Too many arguments"
            
        write(unitReport, *) "Initial condition :"
        
        select case (init)
            case ("rand")
                call randomDeposition(sph_X, sph_rMin)
                write(unitReport, *) "    Random deposition"
            case default
                write(*, *) "Enter the initial condition : "
                write(*, *) "   'rand'."
                stop
        end select
        
    end subroutine initialCondition
    
    !> Primitive cubic configuration
    
    subroutine primitiveCubic(sph_X, sph_rMin)
    
        real(DP), dimension(:, :), intent(inout) :: sph_X
        real(DP), intent(in) :: sph_rMin
    
        integer :: iDir
        integer :: i, j, k, iCol
        integer :: sph_Ncol
        integer, dimension(Dim) :: nCols
        real(DP), dimension(Dim) :: ratio
        real(DP) :: oneThird = 1._DP/3._DP
        
        write(*, *) "Primitive cubic"
        
        ! Proportion according to the direction
        sph_Ncol = size(sph_X, 2)
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
            if ( sph_rMin*real(nCols(iDir), DP) > Lsize(iDir) ) then
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
            sph_X(iDir, :) = sph_X(iDir, :) - 0.5_DP*ratio(iDir) 
            ! just inside
        end do
    
    end subroutine primitiveCubic
    
    !> Random deposition configuration
    
    subroutine randomDeposition(sph_X, sph_rMin)
    
        real(DP), dimension(:, :), intent(inout) :: sph_X
        real(DP), intent(in) :: sph_rMin
    
        integer :: iCol, Ncols, nOK
        integer :: sph_Ncol
        real(DP), dimension(Dim) :: xTest
        real(DP) :: rTest
        
        write(*, *) "Random deposition"
    
        call random_number(sph_X(:, 1))
        sph_X(:, 1) = sph_X(:, 1)*(Lsize(:)-sph_rMin)
        Ncols = 1        
        
        sph_Ncol = size(sph_X, 2)
        do while (Ncols<sph_Ncol)
        
            call random_number(xTest)
            xTest(:) = xTest(:)*(Lsize(:)-sph_rMin)
            
            nOK = 0
            do iCol = 1, Ncols
                rTest = dist(sph_X(:, iCol), xTest(:))
                if (rTest >= sph_rMin) then
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
            sph_X(:, iCol) = sph_X(:, iCol) + inter_radius
        end do
    
    end subroutine randomDeposition
    
    subroutine report(unitReport)
    
    	integer, intent(in) :: unitReport
    
    	write(unitReport, *) "Simulation MC_C :"
        write(unitReport ,*) "    Lsize(:) = ", Lsize(:)
        write(unitReport ,*) "    Vol = ", product(Lsize)
        write(unitReport, *) "    Nstep = ", Nstep
        write(unitReport, *) "    Ntherm = ", Ntherm
        write(unitReport, *) "    Nmove = ", Nmove
    
    end subroutine report
    
end module mod_tools
