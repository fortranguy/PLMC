!> \brief Subroutines for the main program

module module_arguments_monte_carlo

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use module_types_micro, only: Argument_Random, Argument_Initial, Monte_Carlo_Arguments
use module_arguments, only: arg_to_file

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
        write(output_unit, *) "    -i, --initial    CONFIGURATION: choose initial configuration"
        write(output_unit, *) "                     CONFIGURATION = 'r', 'random' (default)"
        write(output_unit, *) "                     CONFIGURATION = 'f', 'files' [dipol_positions] "
        write(output_unit, *) "                     [dipol_orientations] [hardS_positions]"
        write(output_unit, *) "    -r, --random     SEED: choose random number generator seed"
        write(output_unit, *) "                     SEED = 'v', 'variable' (default)"
        write(output_unit, *) "                     SEED = 'f', 'fix'"
        write(output_unit, *) "                     SEED = 'p', 'put' [size] [seed_1]...[seed_n]"
        write(output_unit, *) "                            Warning: compiler-dependent"
        write(output_unit, *)
        write(output_unit, *) "Report bugs to <salomon.chung@u-pe.fr>."
    
    end subroutine write_help

    subroutine read_init_files(i_arg, arg_init)

        integer, intent(inout) :: i_arg
        type(Argument_Initial), intent(inout) :: arg_init

        character(len=4096) :: file
        integer :: i_file, length

        do i_file = 1, 3
            i_arg = i_arg + 1
            call arg_to_file(i_arg, file, length)
            arg_init%files(i_file) = file
            arg_init%length(i_file) = length
        end do
        
    end subroutine read_init_files

    subroutine read_seed_put(i_arg, arg_rand)

        integer, intent(inout) :: i_arg
        type(Argument_Random), intent(inout) :: arg_rand
        
        character(len=4096) :: argument
        integer :: length, status
        integer :: seed_size, arg_rand_size
        integer :: iSeed, arg_rand_i
        
        i_arg = i_arg + 1
        call get_command_argument(i_arg, argument, length, status)
        if (status /= 0) error stop "Error: read_seed_put"
        read(argument(1:length), '(i3)') arg_rand_size ! limits ?

        call random_seed(size = seed_size)
        if (arg_rand_size /= seed_size) error stop "error seed size"

        allocate(arg_rand%seed(seed_size))

        do iSeed = 1, seed_size
            i_arg = i_arg + 1
            call get_command_argument(i_arg, argument, length, status)
            if (status /= 0) error stop "Error: read_seed_put: component"
            read(argument(1:length), '(i11)') arg_rand_i ! limits ?
            arg_rand%seed(iSeed) = arg_rand_i
        end do

    end subroutine read_seed_put

    !> Read arguments

    subroutine read_arguments(args)
    
        type(Monte_Carlo_Arguments), intent(out) :: args

        character(len=4096) :: argument, sub_argument
        integer :: i_arg, length, status
        logical :: rand_redefined, init_redefined

        args%random%choice = 'v'
        args%initial%choice = 'r'

        rand_redefined = .false.
        init_redefined = .false.

        i_arg = 1
        do while(i_arg <= command_argument_count())

            call get_command_argument(i_arg, argument, length, status)
            if (status /= 0) error stop "Error: read_arguments"

            select case (argument)

                case ("-h", "--help")
                    call write_help()
                    stop

                case ("-i", "--initial")
                    if (init_redefined) error stop "Error: initial configuration already defined."
                    i_arg = i_arg + 1
                    call get_command_argument(i_arg, sub_argument, length, status)
                    if (status /= 0) error stop "Enter initial configuration, cf. help."
                    select case (sub_argument)
                        case ("r", "random")
                        case ("f", "files")
                            args%initial%choice = 'f'
                            call read_init_files(i_arg, args%initial)
                        case default
                            call write_help()
                            error stop
                    end select
                    init_redefined = .true.

                case ("-r", "--random")
                    if (rand_redefined) error stop "Error: random seed already defined."
                    i_arg = i_arg + 1
                    call get_command_argument(i_arg, sub_argument, length, status)
                    if (status /= 0) error stop "Enter random seed choice, cf. help."
                    select case (sub_argument)
                        case ("v", "variable")
                        case ("f", "fix")
                            args%random%choice = 'f'
                        case ("p", "put")
                            args%random%choice = 'p'
                            call read_seed_put(i_arg, args%random)
                        case default
                            call write_help()
                            error stop
                    end select
                    rand_redefined = .true.

                case default
                    write(error_unit, *) "Unknown option: '", trim(argument), "'"
                    error stop

            end select

            i_arg = i_arg + 1

        end do

    end subroutine read_arguments

end module module_arguments_monte_carlo
