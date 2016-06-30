module procedures_random_seed_factory

use json_module, only: json_core, json_file, json_value
use procedures_errors, only: error_exit
use procedures_checks, only: check_array_size, check_data_found

implicit none

private
public :: set, write

contains

    subroutine set(input_data, prefix)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: seed_name, data_field
        integer, allocatable :: custom_seed(:)
        logical :: data_found

        data_field = prefix//"name"
        call input_data%get(data_field, seed_name, data_found)
        call check_data_found(data_field, data_found)
        select case (seed_name)
            case ("default")
                call random_seed()
            case ("urandom")
                call set_from_urandom()
            case ("custom")
                data_field = prefix//"seed"
                call input_data%get(data_field, custom_seed, data_found)
                call check_data_found(data_field, data_found)
                call set_from_seed(custom_seed)
            case default
                call error_exit("procedures_random_seed_factory: set: seed_name type unknown. "//&
                    "Choose between: 'default', 'urandom', 'custom'.")
        end select
    end subroutine set

    !> cf. http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
    subroutine set_from_urandom()

        integer, allocatable :: seed(:)
        integer :: seed_size, urandom_unit, i_stat

        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        !Try if the OS provides a random number generator
        open(newunit=urandom_unit, file="/dev/urandom", access="stream", form="unformatted", &
            action="read", status="old", iostat=i_stat)
        if (i_stat /= 0) then
            call error_exit("procedures_random_number: set_from_urandom: "//&
                "I can't open /dev/urandom.")
        end if
        read(urandom_unit) seed
        close(urandom_unit)
        call random_seed(put=seed)
    end subroutine set_from_urandom

    subroutine set_from_seed(custom_seed)
        integer, intent(in) :: custom_seed(:)

        integer :: seed_size

        call random_seed(size=seed_size)
        call check_array_size("procedures_random_seed_factory: set_from_seed", &
            "custom_seed", custom_seed, seed_size)
        call random_seed(put=custom_seed)
    end subroutine set_from_seed

    subroutine write(json, output_data, object_name)
        type(json_core), intent(inout) :: json
        type(json_value), intent(inout), pointer :: output_data
        character(len=*), intent(in) :: object_name

        type(json_value), pointer :: seed_data => null()
        integer, allocatable :: seed(:)
        integer :: seed_size

        call json%create_object(seed_data, object_name)
        call json%add(output_data, seed_data)

        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        call random_seed(get=seed)
        call json%add(seed_data, "seed", seed)
        seed_data => null()
    end subroutine write

end module procedures_random_seed_factory
