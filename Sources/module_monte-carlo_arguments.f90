!> \brief Subroutines for the main program

module module_monteCarlo_arguments

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use module_types, only: argument_seed, argument_initial

implicit none

private
public read_arguments

contains

    subroutine write_help()
    
        write(output_unit, *) "Usage: mc_[ENSEMBLE]_[GEOMETRY] [OPTION]"
        write(output_unit, *) "Run ENSEMBLE Monte-Carlo simulation in GEOMETRY."
        write(output_unit, *)
        write(output_unit, *) "Mandatory arguments to long options are mandatory for short options too."
        write(output_unit, *) "    -h, --help"
        write(output_unit, *) "    -i, --intial CONDITION   CONDITION='r', 'random' (default)"
        write(output_unit, *) "                             CONDITION='f', 'files' [dipol_positions] "
        write(output_unit, *) "                             [dipol_orientations] [hardS_positions]"
        write(output_unit, *) "    -s, --seed STATE         STATE='v', 'variable' (default)"
        write(output_unit, *) "                             STATE='f', 'fix'"
        write(output_unit, *) "                             STATE='p', 'put' [size] [seed_1]...[seed_n]"
        write(output_unit, *) "                                              Warning: compiler-dependent"
        write(output_unit, *)
        write(output_unit, *) "Report bugs to <salomon.chung@u-pe.fr>."
    
    end subroutine write_help

    subroutine read_init_files(iArg, arg_init)

        integer, intent(inout) :: iArg
        type(argument_initial), intent(inout) :: arg_init

        character(len=4096) :: file
        integer :: iFile, length, status
        logical :: exist

        do iFile = 1, 3
            iArg = iArg + 1
            call get_command_argument(iArg, file, length, status)
            if (status /= 0) error stop "no file"
            inquire(file=file(1:length), exist=exist)
            if (.not.exist) then
                write(error_unit, *) "missing file: ", file(1:length)
                error stop
            end if
            arg_init%files(iFile) = file
            arg_init%length(iFile) = length
        end do
        
    end subroutine read_init_files

    subroutine read_seed_put(iArg, arg_seed)

        integer, intent(inout) :: iArg
        type(argument_seed), intent(inout) :: arg_seed
        
        character(len=4096) :: argument
        integer :: length, status
        integer :: seed_size, arg_seed_size
        integer :: iSeed, arg_seed_i
        
        iArg = iArg + 1
        call get_command_argument(iArg, argument, length, status)
        if (status /= 0) error stop "Error: read_seed_put"
        read(argument(1:length), '(i3)') arg_seed_size ! limits ?

        call random_seed(size = seed_size)
        if (arg_seed_size /= seed_size) error stop "error seed size"

        allocate(arg_seed%seed(seed_size))

        do iSeed = 1, seed_size
            iArg = iArg + 1
            call get_command_argument(iArg, argument, length, status)
            if (status /= 0) error stop "Error: read_seed_put: component"
            read(argument(1:length), '(i11)') arg_seed_i ! limits ?
            arg_seed%seed(iSeed) = arg_seed_i
        end do

    end subroutine read_seed_put

    !> Read arguments

    subroutine read_arguments(arg_seed, arg_init)
    
        type(argument_seed), intent(out) :: arg_seed
        type(argument_initial), intent(out) :: arg_init

        character(len=4096) :: argument, sub_argument
        integer :: iArg, length, status
        logical :: seed_redefined, init_redefined

        arg_seed%choice = 'v'
        arg_init%choice = 'r'

        seed_redefined = .false.
        init_redefined = .false.

        iArg = 1
        do while(iArg <= command_argument_count())

            call get_command_argument(iArg, argument, length, status)
            if (status /= 0) error stop "Error: read_arguments"

            select case (argument)

                case ("-h", "--help")
                    call write_help()
                    stop

                case ("-i", "--initial")
                    if (init_redefined) error stop "Error: initial condition already defined."
                    iArg = iArg + 1
                    call get_command_argument(iArg, sub_argument, length, status)
                    if (status /= 0) error stop "Enter initial condition, cf. help."
                    select case (sub_argument)
                        case ("r", "random")
                        case ("f", "files")
                            arg_init%choice = 'f'
                            call read_init_files(iArg, arg_init)
                        case default
                            call write_help()
                            error stop
                    end select
                    init_redefined = .true.

                case ("-s", "--seed")
                    if (seed_redefined) error stop "Error: seed already defined."
                    iArg = iArg + 1
                    call get_command_argument(iArg, sub_argument, length, status)
                    if (status /= 0) error stop "Enter seed choice, cf. help."
                    select case (sub_argument)
                        case ("v", "variable")
                        case ("f", "fix")
                            arg_seed%choice = 'f'
                        case ("p", "put")
                            arg_seed%choice = 'p'
                            call read_seed_put(iArg, arg_seed)
                        case default
                            call write_help()
                            error stop
                    end select
                    seed_redefined = .true.

                case default
                    write(error_unit, *) "Unknown option: '", trim(argument), "'"
                    error stop

            end select

            iArg = iArg + 1

        end do

    end subroutine read_arguments

end module module_monteCarlo_arguments
