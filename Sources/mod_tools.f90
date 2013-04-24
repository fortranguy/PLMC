!> \brief Subroutines as tools for the main program

module mod_tools

use data_constants
use data_particles
use data_mc
use mod_physics
use class_spheres

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
    
    subroutine initialCondition(type1, type2, mix_rMin, unitReport)
    
        class(Spheres), intent(inout) :: type1, type2
        real(DP), intent(in) :: mix_rMin
        integer, intent(in) :: unitReport        

        character(len=20) :: init
        integer :: longueur, statut
        
        call get_command_argument(1, init, longueur, statut)
        if (statut /= 0) stop "error get_command_argument"
        if (command_argument_count() > 1) stop "Too many arguments"
            
        write(unitReport, *) "Initial condition :"
        
        select case (init)
            case ("rand")
                call randomDeposition(type1%X, type1%getRMin(), type2%X, &
                    type2%getRmin(), mix_rMin)
                write(unitReport, *) "    Random deposition"
            case default
                write(*, *) "Enter the initial condition : "
                write(*, *) "   'rand'."
                stop
        end select
        
    end subroutine initialCondition
    
    !> Random deposition configuration
    
    subroutine randomDeposition(type1_X, type1_rMin, type2_X, type2_rMin, &
        mix_rMin)
    
        real(DP), dimension(:, :), intent(inout) :: type1_X, type2_X
        real(DP), intent(in) :: type1_rMin, type2_rMin, mix_rMin
        
        integer :: type1_Ncol, type2_Ncol
        integer :: iCol, iColTest
        real(DP), dimension(Dim) :: xRand
        real(DP) :: rTest
        
        write(*, *) "Random deposition"
        
        ! Type 1
        type1_Ncol = size(type1_X, 2)
        do iCol = 1, type1_Ncol
        
7101        call random_number(xRand)
            type1_X(:, iCol) = xRand*Lsize(:)
            
            do iColTest = 1, iCol-1            
                rTest = dist(type1_X(:, iColTest), type1_X(:, iCol))
                if (rTest < type1_rMin) then
                    goto 7101
                end if            
            end do
        
        end do
        
        ! Type 2
        type2_Ncol = size(type2_X, 2)
        do iCol = 1, type2_Ncol
        
7102        call random_number(xRand)
            type2_X(:, iCol) = xRand*Lsize(:)
            
            do iColTest = 1, type1_Ncol
                rTest = dist(type1_X(:, iColTest), type2_X(:, iCol))
                if (rTest < mix_rMin) then
                    goto 7102
                end if            
            end do
            
            do iColTest = 1, iCol-1            
                rTest = dist(type2_X(:, iColTest), type2_X(:, iCol))
                if (rTest < type2_rMin) then
                    goto 7102
                end if            
            end do
        
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
    
    subroutine results(ePot_mc, ePot_total, ePot_mcSum, duration, unitReport)
    
        real(DP), intent(in) :: ePot_mc, ePot_total, ePot_mcSum
        real(DP), intent(in) :: duration
        integer, intent(in) :: unitReport
        
        write(unitReport, *) "Consistency test:"
        write(unitReport, *) "    ePot_mc = ", ePot_mc
        write(unitReport, *) "    ePot_final = ", ePot_total
        write(unitReport, *) "    relative difference = ", &
            abs(ePot_total-ePot_mc)/ePot_total
            
        write(unitReport, *) "Results :"
        write(unitReport, *) "    average energy = ", &
            ePot_mcSum/real(Nstep, DP)
        write(unitReport, *) "    duration =", duration/60._DP, "min"
    
    end subroutine results
    
end module mod_tools