!> \brief Subroutines as tools for the main program

module mod_tools

use iso_fortran_env
use data_constants
use data_particles
use data_mc
use mod_physics
use class_spheres

implicit none

contains

    !> Random number generator : seed
    
    subroutine init_random_seed(report_unit)
    
        integer, intent(in) :: report_unit
    
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed(:) = clock + 37 * [ (i - 1, i = 1, n) ]
        call random_seed(put = seed)
        
        write(report_unit, *) "RNG :"
        write(report_unit ,*) "    n = ", n
        write(report_unit ,*) "    seed(:) = ", seed(:)

        deallocate(seed)
        
    end subroutine init_random_seed
    
    !> Initial condition
    
    subroutine initialCondition(type1, type2, mix_rMin, report_unit)
    
        class(Spheres), intent(inout) :: type1, type2
        real(DP), intent(in) :: mix_rMin
        integer, intent(in) :: report_unit        

        character(len=20) :: init
        integer :: longueur, statut
        
        call get_command_argument(1, init, longueur, statut)
        if (statut /= 0) stop "error get_command_argument"
        if (command_argument_count() > 1) stop "Too many arguments"
            
        write(report_unit, *) "Initial condition :"
        
        select case (init)
            case ("rand")
                call randomDeposition(type1%X, type1%getRMin(), type2%X, &
                    type2%getRmin(), mix_rMin)
                write(report_unit, *) "    Random deposition"
            case default
                write(output_unit, *) "Enter the initial condition : "
                write(output_unit, *) "   'rand'."
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
        
        write(output_unit, *) "Random deposition"
        
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
    
    !> Report : general
    
    subroutine report(report_unit)
    
        integer, intent(in) :: report_unit

        write(report_unit, *) "Simulation MC_C :"
        write(report_unit ,*) "    Lsize(:) = ", Lsize(:)
        write(report_unit ,*) "    Vol = ", product(Lsize)
        write(report_unit, *) "    Nstep = ", Nstep
        write(report_unit, *) "    Ntherm = ", Ntherm
        write(report_unit, *) "    Nmove = ", Nmove
    
    end subroutine report
    
    !> Results : general
    
    subroutine results(ePot_mc, ePot_total, ePot_mcSum, duration, report_unit)
    
        real(DP), intent(in) :: ePot_mc, ePot_total, ePot_mcSum
        real(DP), intent(in) :: duration
        integer, intent(in) :: report_unit
        
        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    ePot_mc = ", ePot_mc
        write(report_unit, *) "    ePot_final = ", ePot_total
        write(report_unit, *) "    relative difference = ", &
            abs((ePot_total-ePot_mc)/ePot_total)
            
        write(report_unit, *) "Results :"
        write(report_unit, *) "    average energy = ", &
            ePot_mcSum/real(Nstep, DP)
        write(report_unit, *) "    duration =", duration/60._DP, "min"
    
    end subroutine results
    
    !> Results : mix
    
    subroutine mix_results(ePot_mc, ePot_total, ePot_mcSum, report_unit)
    
        real(DP), intent(in) :: ePot_mc, ePot_total, ePot_mcSum
        integer, intent(in) :: report_unit

        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    ePot_mc = ", ePot_mc
        write(report_unit, *) "    ePot_final = ", ePot_total
        write(report_unit, *) "    relative difference = ", &
            abs((ePot_total-ePot_mc)/ePot_total)
    
        write(report_unit, *) "Results :"        
        write(report_unit, *) "    average energy = ", ePot_mcSum/real(Nstep, DP)
    
    end subroutine mix_results
    
end module mod_tools
