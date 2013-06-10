!> \brief Subroutines as tools for the main program

module mod_tools

use, intrinsic :: iso_fortran_env
use data_constants
use data_particles
use data_mc
use mod_physics
use class_spheres
use class_dipolarSpheres
use class_hardSpheres

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
    
        class(DipolarSpheres), intent(inout) :: type1
        class(Spheres), intent(inout) :: type2
        real(DP), intent(in) :: mix_rMin
        integer, intent(in) :: report_unit        

        character(len=4096) :: init
        integer :: length, status
        
        write(report_unit, *) "Initial condition :"
        
        select case (command_argument_count())
        
            case (1)
            
                call get_command_argument(1, init, length, status)
                if (status /= 0) stop "error get_command_argument"
                
                select case (init)
                    case ("rand") 
                        call randomDepositions(type1%positions, type1%getRMin(), type2%positions, &
                                               type2%getRmin(), mix_rMin)
                        call randomorientations(type1%orientations)
                        write(output_unit, *) "Random depositions + random orientations"
                        write(report_unit, *) "    Random depositions + random orientations"
                    case default
                        write(error_unit, *) "Enter the initial condition : "
                        write(error_unit, *) &
                            "   'rand' or '[type1_positions] [type1_orientations] [type2_positions]'."
                        stop
                end select
                
            case (3)
            
                write(output_unit, *) "Old configuration"
                write(report_unit, *) "    Old configuration"
                
                call oldConfiguration(1, type1%getName()//"_positions", type1%positions, &
                                      dot_product(Lsize, Lsize))
                call oldConfiguration(2, type1%getName()//"_orientations", type1%orientations, 1._DP)
                call oldConfiguration(3, type2%getName()//"_positions", type2%positions, &
                                      dot_product(Lsize, Lsize))
            
            case default
                write(error_unit, *) "Enter the initial condition : "
                write(error_unit, *) &
                    "   'rand' or '[type1_positions] [type1_orientations] [type2_positions]'."
                stop
                
        end select
        
    end subroutine initialCondition
    
    !> Random depositions configuration
    
    subroutine randomDepositions(type1_positions, type1_rMin, type2_positions, type2_rMin, mix_rMin)
    
        real(DP), dimension(:, :), intent(out) :: type1_positions, type2_positions
        real(DP), intent(in) :: type1_rMin, type2_rMin, mix_rMin
        
        integer :: type1_Ncol, type2_Ncol
        integer :: iCol, iColTest
        real(DP), dimension(Dim) :: xRand
        real(DP) :: rTest
        
        ! Type 1
        type1_Ncol = size(type1_positions, 2)
        do iCol = 1, type1_Ncol
        
7101        call random_number(xRand)
            type1_positions(:, iCol) = xRand*Lsize(:)
            
            do iColTest = 1, iCol-1            
                rTest = dist(type1_positions(:, iColTest), type1_positions(:, iCol))
                if (rTest < type1_rMin) then
                    goto 7101
                end if            
            end do
        
        end do
        
        ! Type 2
        type2_Ncol = size(type2_positions, 2)
        do iCol = 1, type2_Ncol
        
7102        call random_number(xRand)
            type2_positions(:, iCol) = xRand*Lsize(:)
            
            do iColTest = 1, type1_Ncol
                rTest = dist(type1_positions(:, iColTest), type2_positions(:, iCol))
                if (rTest < mix_rMin) then
                    goto 7102
                end if            
            end do
            
            do iColTest = 1, iCol-1            
                rTest = dist(type2_positions(:, iColTest), type2_positions(:, iCol))
                if (rTest < type2_rMin) then
                    goto 7102
                end if            
            end do
        
        end do
    
    end subroutine randomDepositions
    
    !> Uniform (gaussian) orientations
    
    subroutine randomorientations(type_orientations)
    
        real(DP), dimension(:, :), intent(out) :: type_orientations
        
        integer :: type_Ncol
        integer :: iCol
        
        type_Ncol = size(type_orientations, 2)
        
        do iCol = 1, type_Ncol
            type_orientations(:, iCol) = random_surface()
        end do
    
    end subroutine randomorientations
    
    !> From an old configuration
    
    subroutine oldConfiguration(iFile, type_name, type_vec, normSqrMax)
    
        integer, intent(in) :: iFile
        character(len=*), intent(in) :: type_name
        real(DP), dimension(:, :), intent(out) :: type_vec
        real(DP), intent(in) :: normSqrMax
    
        character(len=20) :: file
        integer :: length
        integer :: fileStat, readStat
        integer :: file_unit

        integer :: type_Ncol
        integer :: iCol
        real(DP), dimension(Dim) :: vecDummy
        real(DP) :: normSqr
        
        call get_command_argument(iFile, file, length, fileStat)
        if (fileStat /= 0) stop "error get_command_argument"
        write(output_unit, *) type_name, " <- ", file(1:length)
        open(newunit=file_unit, recl=4096, file=file(1:length), status='old', action='read')
        
        type_Ncol = size(type_vec, 2)
        
        iCol = 0
        do
            read(file_unit, fmt=*, iostat=readStat) vecDummy(:)
            if (readStat == iostat_end) exit
            iCol = iCol + 1
        end do
        
        if (iCol == type_Ncol) then
            rewind(file_unit)
            do iCol = 1, type_Ncol
                read(file_unit, *) type_vec(:, iCol)
                normSqr = dot_product(type_vec(:, iCol), type_vec(:, iCol))
                if (normSqr > normSqrMax+io_tiny) then
                    write(error_unit, *) "Norm error : ", file(1:length)
                    write(error_unit, *) "Vec ", type_vec(:, iCol)
                    write(error_unit, *) "NormSqr", normSqr
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
        write(report_unit ,*) "    Volume = ", Volume
        write(report_unit ,*) "    kMax(:) = ", kMax(:)
        write(report_unit ,*) "    NwaveVectors =", (2*kMax(1)+1) * (2*kMax(2)+1) * (2*kMax(3)+1)
        write(report_unit ,*) "    Ncol = ", Ncol
        write(report_unit ,*) "    Tstar = ", Tstar
        
        write(report_unit, *) "    Nstep = ", Nstep
        write(report_unit, *) "    Ntherm = ", Ntherm
        write(report_unit, *) "    Nmove = ", Nmove
        write(report_unit, *) "    Nrotate = ", Nrotate
    
    end subroutine report
    
    ! Total & Mix : consistency test
    
    subroutine consistTest(Epot, Epot_conf, report_unit)
    
        real(DP), intent(in) :: Epot, Epot_conf
        integer, intent(in) :: report_unit
        
        real(DP) :: difference
    
        difference = abs((Epot_conf-Epot)/Epot_conf)
        
        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    Epot = ", Epot
        write(report_unit, *) "    Epot_conf = ", Epot_conf
        write(report_unit, *) "    relative difference = ", difference
        
        if (difference > consist_tiny) then
            write(report_unit, *) "    WARNING !"
        else
            write(report_unit, *) "    OK !"
        end if
    
    end subroutine consistTest
    
    !> Total : Results    
    
    subroutine printResults(EpotSum, duration, report_unit)
    
        real(DP), intent(in) :: EpotSum
        real(DP), intent(in) :: duration
        integer, intent(in) :: report_unit
            
        write(report_unit, *) "Results :"
        write(report_unit, *) "    average energy = ", EpotSum/real(Nstep, DP)
        write(report_unit, *) "    average energy per particule = ", &
                                   EpotSum/real(Nstep, DP)/real(Ncol, DP)
        write(report_unit, *) "    duration =", duration/60._DP, "min"
    
    end subroutine printResults
    
    !> Mix : Results
    
    subroutine mix_printResults(EpotSum, report_unit)
    
        real(DP), intent(in) :: EpotSum
        integer, intent(in) :: report_unit
    
        write(report_unit, *) "Results :"        
        write(report_unit, *) "    average energy = ", EpotSum/real(Nstep, DP)
    
    end subroutine mix_printResults
    
end module mod_tools
