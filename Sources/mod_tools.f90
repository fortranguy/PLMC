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
    
    subroutine initRandomSeed(report_unit)
    
        integer, intent(in) :: report_unit
    
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed(:) = clock + 37 * [ (i - 1, i = 1, n) ]
        call random_seed(put = seed)
        
        write(report_unit, *) "Random number generator :"
        write(report_unit ,*) "    n = ", n
        write(report_unit ,*) "    seed(:) = ", seed(:)

        deallocate(seed)
        
    end subroutine initRandomSeed
    
    !> Initial condition
    
    subroutine initialCondition(type1, type2, mix_rMin, report_unit)
    
        class(Spheres), intent(inout) :: type1, type2
        real(DP), intent(in) :: mix_rMin
        integer, intent(in) :: report_unit        

        character(len=20) :: init
        integer :: length, status
        
        write(report_unit, *) "Initial condition :"
        
        select case (command_argument_count())
        
            case (1)
            
                call get_command_argument(1, init, length, status)
                if (status /= 0) stop "error get_command_argument"
                
                select case (init)
                    case ("rand") 
                        call randomDeposition(type1%X, type1%getRMin(), type2%X, type2%getRmin(), mix_rMin)
                        write(output_unit, *) "Random deposition"
                        write(report_unit, *) "    Random deposition"
                    case default
                        write(error_unit, *) "Enter the initial condition : "
                        write(error_unit, *) "   'rand' or '[file1] [file2]'."
                        stop
                end select
                
            case (2)
            
                write(output_unit, *) "Old configuration"
                write(report_unit, *) "    Old configuration"
                
                call oldConfiguration(1, type1%getName(), type1%X)
                call oldConfiguration(2, type2%getName(), type2%X)
            
            case default
                write(error_unit, *) "Enter the initial condition : "
                write(error_unit, *) "   'rand' or '[file1] [file2]'."
                stop
                
        end select
        
    end subroutine initialCondition
    
    !> Random deposition configuration
    
    subroutine randomDeposition(type1_X, type1_rMin, type2_X, type2_rMin, mix_rMin)
    
        real(DP), dimension(:, :), intent(inout) :: type1_X, type2_X
        real(DP), intent(in) :: type1_rMin, type2_rMin, mix_rMin
        
        integer :: type1_Ncol, type2_Ncol
        integer :: iCol, iColTest
        real(DP), dimension(Dim) :: xRand
        real(DP) :: rTest
        
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
    
    !> From an old configuration
    
    subroutine oldConfiguration(iType, type_name, type_X)
    
        integer, intent(in) :: iType
        character(len=*), intent(in) :: type_name
        real(DP), dimension(:, :), intent(out) :: type_X
    
        character(len=20) :: file
        integer :: length
        integer :: fileStat, readStat
        integer :: file_unit

        integer :: type_Ncol
        integer :: iCol
        real(DP), dimension(Dim) :: xDummy
        real(DP) :: normSqr, normSqrMax
        
        call get_command_argument(iType, file, length, fileStat)
        if (fileStat /= 0) stop "error get_command_argument"
        write(output_unit, *) type_name, " <- ", file(1:length)
        open(newunit=file_unit, recl=4096, file=file(1:length), status='old', action='read')
        
        type_Ncol = size(type_X, 2)
        normSqrMax = dot_product(Lsize, Lsize)
        
        iCol = 0
        do
            read(file_unit, fmt=*, iostat=readStat) xDummy(:)
            if (readStat == iostat_end) exit
            iCol = iCol + 1
        end do
        
        if (iCol == type_Ncol) then
            rewind(file_unit)
            do iCol = 1, type_Ncol
                read(file_unit, *) type_X(:, iCol)
                normSqr = dot_product(type_X(:, iCol), type_X(:, iCol))
                if (normSqr > normSqrMax) then
                    write(error_unit, *) "Size error : ", file(1:length)
                    write(error_unit, *) "xCol ", type_X(:, iCol)
                    write(error_unit, *) "vs"
                    write(error_unit, *) "Lsize", Lsize(:)
                    stop
                end if
            end do
        else
            write(error_unit, *) "Error reading : ", file(1:length)
            write(error_unit, *) "iCol", iCol, " /= ", "type_Ncol", type_Ncol
            stop
        end if        
        
        close(file_unit)
        
    end subroutine oldConfiguration
    
    !> Total : report
    
    subroutine report(report_unit)
    
        integer, intent(in) :: report_unit

        write(report_unit, *) "Data :"
        
        write(report_unit ,*) "    Lsize(:) = ", Lsize(:)
        write(report_unit ,*) "    Vol = ", product(Lsize)
        write(report_unit ,*) "    Ncol = ", Ncol
        write(report_unit ,*) "    Tstar = ", Tstar
        
        write(report_unit, *) "    Nstep = ", Nstep
        write(report_unit, *) "    Ntherm = ", Ntherm
        write(report_unit, *) "    Nmove = ", Nmove
    
    end subroutine report
    
    ! Total : consistency test
    
    subroutine consistTest(Epot, Epot_conf, report_unit)
    
        real(DP), intent(in) :: Epot, Epot_conf
        integer, intent(in) :: report_unit
    
        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    Epot = ", Epot
        write(report_unit, *) "    Epot_conf = ", Epot_conf
        write(report_unit, *) "    relative difference = ", abs((Epot_conf-Epot)/Epot_conf)
    
    end subroutine consistTest
    
    !> Total : Results    
    
    subroutine results(EpotSum, duration, report_unit)
    
        real(DP), intent(in) :: EpotSum
        real(DP), intent(in) :: duration
        integer, intent(in) :: report_unit
            
        write(report_unit, *) "Results :"
        write(report_unit, *) "    average energy = ", EpotSum/real(Nstep, DP)
        write(report_unit, *) "    duration =", duration/60._DP, "min"
    
    end subroutine results
    
    ! Mix : consistency test
    
    subroutine mix_consistTest(Epot, Epot_conf, report_unit)
    
        real(DP), intent(in) :: Epot, Epot_conf
        integer, intent(in) :: report_unit
    
        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    Epot = ", Epot
        write(report_unit, *) "    Epot_conf = ", Epot_conf
        write(report_unit, *) "    relative difference = ", abs((Epot_conf-Epot)/Epot_conf)
    
    end subroutine mix_consistTest
    
    !> Mix : Results
    
    subroutine mix_results(EpotSum, report_unit)
    
        real(DP), intent(in) :: EpotSum
        integer, intent(in) :: report_unit
    
        write(report_unit, *) "Results :"        
        write(report_unit, *) "    average energy = ", EpotSum/real(Nstep, DP)
    
    end subroutine mix_results
    
end module mod_tools
