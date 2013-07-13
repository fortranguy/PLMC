!> \brief Subroutines as tools for the main program

module module_tools

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, iostat_end
use data_precisions, only : DP, io_tiny, consist_tiny
use data_box, only : Ndim, Lsize, Volume, Kmax
use data_particles, only : Ncol
use data_monteCarlo, only : Temperature, Nstep, decorrelFactor, Ntherm, Nmove, Nrotate
use module_physics, only : dist_PBC, random_surface
use class_spheres
use class_dipolarSpheres
use class_hardSpheres

implicit none
private
public initRandomSeed, initialCondition, report, consistTest, printResults, mix_printResults

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
                        call randomDepositions(type1, type2, mix_rMin)
                        call randomOrientations(type1%orientations, type1%getNcol())
                        write(output_unit, *) "Random depositions + random orientations"
                        write(report_unit, *) "    Random depositions + random orientations"
                    case default
                        write(error_unit, *) "'", trim(init), "'", " isn't a known argument."
                        write(error_unit, *) "Do you mean 'rand' ?"
                        stop
                end select
                
            case (3)
            
                write(output_unit, *) "Old configuration"
                write(report_unit, *) "    Old configuration"
                
                call oldConfiguration(1, type1%getName()//"_positions", type1%getNcol(), &
                                      type1%positions, norm2(Lsize))
                call oldConfiguration(2, type1%getName()//"_orientations", type1%getNcol(), &
                                      type1%orientations, 1._DP)
                ! Warning : be careful with the unit !
                call oldConfiguration(3, type2%getName()//"_positions", type2%getNcol(),  &
                                      type2%positions, norm2(Lsize))
            
            case default
                write(error_unit, *) "Enter the initial condition : "
                write(error_unit, *) &
                    "   'rand' or '[type1_positions] [type1_orientations] [type2_positions]'."
                stop
                
        end select
        
    end subroutine initialCondition
    
    !> Random depositions configuration
    
    subroutine randomDepositions(type1, type2, mix_rMin)

        class(DipolarSpheres), intent(inout) :: type1
        class(Spheres), intent(inout) :: type2
        real(DP), intent(in) :: mix_rMin

        integer :: iCol, iColTest
        real(DP), dimension(Ndim) :: xRand
        real(DP) :: rTest
        
        ! Type 1
        do iCol = 1, type1%getNcol()
        
7101        continue
            call random_number(xRand)
            type1%positions(:, iCol) = xRand*Lsize(:)
            
            do iColTest = 1, iCol-1            
                rTest = dist_PBC(type1%positions(:, iColTest), type1%positions(:, iCol))
                if (rTest < type1%getRmin()) then
                    goto 7101
                end if            
            end do
        
        end do
        
        ! Type 2
        do iCol = 1, type2%getNcol()
        
7102        continue
            call random_number(xRand)
            type2%positions(:, iCol) = xRand*Lsize(:)
            
            do iColTest = 1, type1%getNcol()
                rTest = dist_PBC(type1%positions(:, iColTest), type2%positions(:, iCol))
                if (rTest < mix_rMin) then
                    goto 7102
                end if            
            end do
            
            do iColTest = 1, iCol-1            
                rTest = dist_PBC(type2%positions(:, iColTest), type2%positions(:, iCol))
                if (rTest < type2%getRmin()) then
                    goto 7102
                end if            
            end do
        
        end do
    
    end subroutine randomDepositions
    
    !> Uniform (gaussian) orientations
    
    subroutine randomOrientations(orientations, Ncol)
    
        real(DP), dimension(:, :), intent(out) :: orientations
        integer, intent(in) :: Ncol
        
        integer :: iCol
        
        do iCol = 1, Ncol
            orientations(:, iCol) = random_surface()
        end do
    
    end subroutine randomOrientations
    
    !> From an old configuration
    
    subroutine oldConfiguration(iFile, type_name, type_Ncol, type_coords, normMax)
    
        integer, intent(in) :: iFile
        character(len=*), intent(in) :: type_name
        integer, intent(in) :: type_Ncol
        real(DP), dimension(:, :), intent(out) :: type_coords
        real(DP), intent(in) :: normMax
    
        character(len=20) :: file
        integer :: length
        integer :: fileStat, readStat
        integer :: file_unit

        integer :: iCol
        real(DP), dimension(Ndim) :: vecDummy
        
        call get_command_argument(iFile, file, length, fileStat)
        if (fileStat /= 0) stop "error get_command_argument"
        write(output_unit, *) type_name, " <- ", file(1:length)
        open(newunit=file_unit, recl=4096, file=file(1:length), status='old', action='read')
        
        iCol = 0
        do
            read(file_unit, fmt=*, iostat=readStat) vecDummy(:)
            if (readStat == iostat_end) exit
            iCol = iCol + 1
        end do
        
        if (iCol == type_Ncol) then
            rewind(file_unit)
            do iCol = 1, type_Ncol
                read(file_unit, *) type_coords(:, iCol)
                if (norm2(type_coords(:, iCol)) > normMax+io_tiny) then
                    write(error_unit, *) "Norm error : ", file(1:length)
                    write(error_unit, *) "Coordinates ", type_coords(:, iCol)
                    write(error_unit, *) "Norm", norm2(type_coords(:, iCol))
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
        write(report_unit ,*) "    Kmax(:) = ", Kmax(:)
        write(report_unit ,*) "    NwaveVectors =", (2*Kmax(1)+1) * (2*Kmax(2)+1) * (2*Kmax(3)+1)
        write(report_unit ,*) "    Ncol = ", Ncol
        write(report_unit ,*) "    Temperature = ", Temperature
        
        write(report_unit, *) "    Nstep = ", Nstep
        write(report_unit, *) "    Ntherm = ", Ntherm
        write(report_unit, *) "    decorrelFactor = ", decorrelFactor
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
    
end module module_tools