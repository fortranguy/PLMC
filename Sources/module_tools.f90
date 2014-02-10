!> \brief Subroutines as tools for the main program

module module_tools

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, iostat_end
use data_precisions, only : DP, real_zero, io_tiny, consist_tiny
use data_constants, only : PI, sigma3d
use data_box, only : Ndim, Lsize, Kmax
use data_monteCarlo, only : Temperature, Nadapt, Nstep, decorrelFactor, Nthermal
use data_potential, only : print_potential
use module_physics, only : dist_PBC, random_surface
use class_hardSpheres
use class_interactingSpheres
use class_dipolarSpheres
use class_mixingPotential
use class_observables
use class_units

implicit none
private
public read_arguments, open_units, mix_open_units, init_randomSeed, set_initialCondition, &
       print_report, mix_init, mix_final, init, final, adapt_move, adapt_rotate, test_consist, &
       print_results, mix_print_results

contains

    subroutine print_help()
    
        write(output_unit, *) "Usage: mc_[ENSEMBLE]_[GEOMETRY] [OPTION]"
        write(output_unit, *) "Run ENSEMBLE Monte-Carlo simulation in GEOMETRY."
        write(output_unit, *)        
        write(output_unit, *) "Mandatory arguments to long options are mandatory for short options too."
        write(output_unit, *) "    -h, --help"
        write(output_unit, *) "    -i, --intial CONDITION   CONDITION='r', 'random' : random desposition"
        write(output_unit, *) "                             CONDITION='f', 'files' [dipol_positions] "
        write(output_unit, *) "                             [dipol_orientations] [hardS_positions]"
        write(output_unit, *) "    -s, --fix-seed           Fix the seed to the default value."
        write(output_unit, *)
        write(output_unit, *) "Report bugs to <salomon.chung@u-pe.fr>."
    
    end subroutine print_help

    !> Read arguments

    subroutine read_arguments(variable_seed)
    
        logical, intent(out) :: variable_seed

        character(len=4096) :: argument
        integer :: iArg, length, status
        
        logical :: rand_initial
        
        rand_initial = .true.
        variable_seed = .true.

        do iArg = 1, command_argument_count()

            call get_command_argument(iArg, argument, length, status)
            if (status /= 0) stop "error get_command_argument"

            select case (argument)

                case ("-h", "--help")
                    call print_help()

                case ("-i", "--initial")
                    call get_command_argument(iArg+1, argument, length, status)
                    if (status /= 0) stop "Enter initial condition, cf. help."
                    select case (argument)
                        case ("r", "random")
                            write(*, *) "random"
                        case ("f", "files")
                            write(*, *) "files"
                        case default
                            call print_help()
                            stop
                    end select

                case ("-s", "--fix-seed")
                    variable_seed = .false.
                    write(*, *) "variable_seed = .false."

            end select

        end do

        stop "end of read_arguments"

    end subroutine read_arguments

    !> Total : open units
    
    subroutine open_units(report_unit, obsThermal_unit, obsEquilib_unit)
    
        integer, intent(out) :: report_unit, obsThermal_unit, obsEquilib_unit
    
        open(newunit=report_unit, recl=4096, file="report.txt", status='new', action='write')
        open(newunit=obsThermal_unit, recl=4096, file="obsThermal.out", status='new', &
             action='write')
        open(newunit=obsEquilib_unit, recl=4096, file="obsEquilib.out", status='new', &
             action='write')
        write(obsEquilib_unit, *) "#", 1 ! 1 observable : energy
        
    end subroutine open_units

    !> Mix : open units
    
    subroutine mix_open_units(mix_report_unit, mix_Epot_unit, mix_obsThermal_unit, &
                                  mix_obsEquilib_unit)
                                  
        integer, intent(out) :: mix_report_unit, mix_Epot_unit, mix_obsThermal_unit, &
                                  mix_obsEquilib_unit
    
        open(newunit=mix_report_unit, recl=4096, file="mix_report.txt", status='new', action='write')
        open(newunit=mix_Epot_unit, recl=4096, file="mix_Epot.tmp", status='new', action='write')
        open(newunit=mix_obsThermal_unit, recl=4096, file="mix_obsThermal.out", &
             status='new', action='write')
        open(newunit=mix_obsEquilib_unit, recl=4096, file="mix_obsEquilib.out", status='new', &
             action='write')
        write(mix_obsEquilib_unit, *) "#", 1 ! 1 observable : energy
        
     end subroutine mix_open_units

    !> Random number generator : seed
    
    subroutine init_randomSeed(variable_seed, report_unit)
    
        logical, intent(in) :: variable_seed
        integer, intent(in) :: report_unit
    
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed
        
        if (variable_seed) then

            call random_seed(size = n)
            allocate(seed(n))

            call system_clock(count=clock)
            
            seed(:) = clock + 37 * [ (i - 1, i = 1, n) ]
            call random_seed(put = seed)
            
            write(report_unit, *) "Random number generator :"
            write(report_unit ,*) "    n = ", n
            write(report_unit ,*) "    seed(:) = ", seed(:)

            deallocate(seed)
            
        end if
        
    end subroutine init_randomSeed
    
    !> Initial condition
    
    subroutine set_initialCondition(dipolar, spherical, mix_sigma, report_unit)
    
        class(DipolarSpheres), intent(inout) :: dipolar
        class(HardSpheres), intent(inout) :: spherical
        real(DP), intent(in) :: mix_sigma
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
                        call randomDepositions(dipolar, spherical, mix_sigma)
                        call randomOrientations(dipolar%orientations, dipolar%get_Ncol())
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
                
                call oldConfiguration(1, dipolar%get_name()//"_positions", dipolar%get_Ncol(), &
                                      dipolar%positions, norm2(Lsize))
                call oldConfiguration(2, dipolar%get_name()//"_orientations", dipolar%get_Ncol(), &
                                      dipolar%orientations, 1._DP)
                call oldConfiguration(3, spherical%get_name()//"_positions", spherical%get_Ncol(), &
                                      spherical%positions, norm2(Lsize))
            
            case default
                write(error_unit, *) "Enter the initial condition : "
                write(error_unit, *) &
                    "   'rand' or '[dipolar_positions] [dipolar_orientations] [spherical_positions]'."
                stop
                
        end select
        
    end subroutine set_initialCondition
    
    !> Random depositions configuration
    
    subroutine randomDepositions(type1, type2, mix_sigma)

        class(HardSpheres), intent(inout) :: type1, type2
        real(DP), intent(in) :: mix_sigma

        integer :: iCol, iColTest
        real(DP), dimension(Ndim) :: xRand
        real(DP) :: rTest
        
        ! Type 1
        do iCol = 1, type1%get_Ncol()
        
7101        continue
            call random_number(xRand)
            type1%positions(:, iCol) = xRand*Lsize(:)
            
            do iColTest = 1, iCol-1
                rTest = dist_PBC(type1%positions(:, iColTest), type1%positions(:, iCol))
                if (rTest < type1%get_sigma()) then
                    goto 7101
                end if
            end do
        
        end do
        
        ! Type 2
        do iCol = 1, type2%get_Ncol()
        
7102        continue
            call random_number(xRand)
            type2%positions(:, iCol) = xRand*Lsize(:)
            
            do iColTest = 1, type1%get_Ncol()
                rTest = dist_PBC(type1%positions(:, iColTest), type2%positions(:, iCol))
                if (rTest < mix_sigma) then
                    goto 7102
                end if
            end do
            
            do iColTest = 1, iCol-1
                rTest = dist_PBC(type2%positions(:, iColTest), type2%positions(:, iCol))
                if (rTest < type2%get_sigma()) then
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
    
    !> Total : print_report
    
    subroutine print_report(Ncol, Nmove, Nrotate, report_unit)
    
        integer, intent(in) :: Ncol, Nmove, Nrotate
        integer, intent(in) :: report_unit

        write(report_unit, *) "Data :"
        
        write(report_unit ,*) "    Precision = ", DP
        write(report_unit ,*) "    Real zero = ", real_zero
        write(report_unit ,*) "    I/O tiny = ", io_tiny
        write(report_unit ,*) "    Energy consistency tiny = ", consist_tiny
        
        write(report_unit ,*) "    Pi = ", PI
        write(report_unit ,*) "    Sigma3d = ", sigma3d
        
        write(report_unit ,*) "    Lsize(:) = ", Lsize(:)
        write(report_unit ,*) "    Volume = ", product(Lsize)
        write(report_unit ,*) "    Kmax(:) = ", Kmax(:)
        write(report_unit ,*) "    NwaveVectors =", (2*Kmax(1)+1) * (2*Kmax(2)+1) * (2*Kmax(3)+1)
        write(report_unit ,*) "    Ncol = ", Ncol
        write(report_unit ,*) "    Temperature = ", Temperature
        
        write(report_unit, *) "    Nstep = ", Nstep
        write(report_unit, *) "    Nthermal = ", Nthermal
        write(report_unit, *) "    decorrelFactor = ", decorrelFactor
        write(report_unit, *) "    Nmove = ", Nmove
        write(report_unit, *) "    Nrotate = ", Nrotate
    
    end subroutine print_report
    
    !> Mix initialisation
    
    subroutine mix_init(mix, type1, type2, mix_Epot_unit, mix_Epot)
    
        class(MixingPotential), intent(inout) :: mix
        class(HardSpheres), intent(in) :: type1, type2
        integer, intent(in) :: mix_Epot_unit
        real(DP), intent(out) :: mix_Epot
    
        call mix%test_overlap(type1, type2)
        call mix%Epot_init()
        if (print_potential) then
            call mix%Epot_print(mix_Epot_unit)
        end if
        call mix%set_cell_size()
        mix_Epot = mix%Epot_conf(type1, type2)
    
    end subroutine mix_init
    
    !> Mix finalization
    
    subroutine mix_final(mix, type1, type2, mix_report_unit, mix_Epot, mix_EpotSum, mix_Epot_conf)
    
        class(MixingPotential), intent(inout) :: mix
        class(HardSpheres), intent(in) :: type1, type2
        integer, intent(in) :: mix_report_unit
        real(DP), intent(in) :: mix_Epot, mix_EpotSum
        real(DP), intent(out) :: mix_Epot_conf
        
        call mix%test_overlap(type1, type2)
        call mix%Epot_init()
        mix_Epot_conf = mix%Epot_conf(type1, type2)
        call test_consist(mix_Epot, mix_Epot_conf, mix_report_unit)
        call mix_print_results(mix_EpotSum, mix_report_unit)
    
    end subroutine mix_final
    
    !> Spheres initialisations
    
    subroutine init(this, other, mix, this_units, this_Epot)
    
        class(HardSpheres), intent(inout) :: this
        class(HardSpheres), intent(in) :: other
        class(MixingPotential), intent(in) :: mix
        class(Units), intent(in) :: this_units
        real(DP), intent(inout) :: this_Epot
        
        call this%test_overlap()
        call this%snap_data(this_units%snap_positions)
        call this%snap_positions(0, this_units%snapIni_positions)
        call this%Epot_init()
        
        if (print_potential) then
            call this%Epot_print(this_units%Epot)
        end if
        select type (this)
            type is (DipolarSpheres)            
                select type (this_units)
                    type is (MoreUnits)
                        call this%snap_data(this_units%snap_orientations)
                        call this%snap_orientations(0, this_units%snapIni_orientations)
                        if (print_potential) then
                            call this%Epot_real_print(this_units%Epot_real)
                        end if
                        call this%Epot_reci_count_waveVectors(this_units%waveVectors)
                end select
        end select        
        this_Epot = this%Epot_conf()        
        
        call this%construct_cells(other, mix%get_cell_size(), mix%get_rCut())
        call this%print_report(this_units%report)
    
    end subroutine init
    
    !> Spheres finalizations
    
    subroutine final(this, this_units, this_obs)
    
        class(HardSpheres), intent(inout) :: this
        class(Units), intent(in) :: this_units
        class(Observables), intent(in) :: this_obs
        
        call this%test_overlap()
        call this%Epot_init()
        call test_consist(this_obs%Epot, this%Epot_conf(), this_units%report)
        call this%snap_positions(0, this_units%snapFin_positions)
        call this_obs%print_results(this_units%report)
        
        select type (this)
            type is (DipolarSpheres)            
                select type (this_units)
                    type is (MoreUnits)
                        call this%snap_orientations(0, this_units%snapFin_orientations)
                end select
        end select
    
    end subroutine final
    
    !> Change : average & adaptation
    
    subroutine adapt_move(this, iStep, obs, move_unit)
    
        class(HardSpheres), intent(inout) :: this
        integer, intent(in) :: iStep
        class(Observables), intent(inout) :: obs
        integer, intent(in) :: move_unit
    
        obs%move_rejectAvg = obs%move_rejectAdapt / real(Nadapt-1, DP)
        obs%move_rejectAdapt = 0._DP
        call this%adapt_move_delta(obs%move_rejectAvg)
        write(move_unit, *) iStep, this%get_move_delta(), obs%move_rejectAvg
    
    end subroutine adapt_move
    
    subroutine adapt_rotate(this, iStep, obs, rotate_unit)
    
        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: iStep
        class(MoreObservables), intent(inout) :: obs
        integer, intent(in) :: rotate_unit
        
        obs%rotate_rejectAvg = obs%rotate_rejectAdapt / real(Nadapt-1, DP)
        obs%rotate_rejectAdapt = 0._DP
        call this%adapt_rotate_delta(obs%rotate_rejectAvg)
        write(rotate_unit, *) iStep, this%get_rotate_delta(), obs%rotate_rejectAvg
        
    end subroutine adapt_rotate
    
    !> Total & Mix : consistency test
    
    subroutine test_consist(Epot, Epot_conf, report_unit)
    
        real(DP), intent(in) :: Epot, Epot_conf
        integer, intent(in) :: report_unit
        
        real(DP) :: difference
        
        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    Epot = ", Epot
        write(report_unit, *) "    Epot_conf = ", Epot_conf
        
        if (abs(Epot_conf) < real_zero) then
            difference = abs(Epot_conf-Epot)
            write(report_unit, *) "    absolute difference = ", difference
        else        
            difference = abs((Epot_conf-Epot)/Epot_conf)
            write(report_unit, *) "    relative difference = ", difference        
        end if
        
        if (difference > consist_tiny) then ! not sufficient for HS ?
            write(report_unit, *) "    WARNING !"
        else
            write(report_unit, *) "    OK !"
        end if
    
    end subroutine test_consist
    
    !> Total : Results
    
    subroutine print_results(Ncol, EpotSum, duration, report_unit)
    
        integer, intent(in) :: Ncol
        real(DP), intent(in) :: EpotSum
        real(DP), intent(in) :: duration
        integer, intent(in) :: report_unit
            
        write(report_unit, *) "Results :"
        write(report_unit, *) "    average energy = ", EpotSum/real(Nstep, DP)
        write(report_unit, *) "    average energy per particule = ", &
                                   EpotSum/real(Nstep, DP)/real(Ncol, DP)
        write(report_unit, *) "    duration =", duration/60._DP, "min"
    
    end subroutine print_results
    
    !> Mix : Results
    
    subroutine mix_print_results(EpotSum, report_unit)
    
        real(DP), intent(in) :: EpotSum
        integer, intent(in) :: report_unit
    
        write(report_unit, *) "Results :"
        write(report_unit, *) "    average energy = ", EpotSum/real(Nstep, DP)
    
    end subroutine mix_print_results
    
end module module_tools
