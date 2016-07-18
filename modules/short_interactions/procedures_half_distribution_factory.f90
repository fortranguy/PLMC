module procedures_half_distribution_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_half_distribution, only: Abstract_Half_Distribution, Rectangular_Half_Distribution, &
    Null_Half_Distribution

implicit none

private
public :: create, destroy

contains

    subroutine create(half_distribution, needed, exploring_data, prefix)
        class(Abstract_Half_Distribution), allocatable, intent(out) :: half_distribution
        logical, intent(in) :: needed
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: distribution_name, data_field
        real(DP) :: delta_distance
        logical :: data_found

        if (needed) then
            data_field = prefix//"Half Distribution.name"
            call exploring_data%get(data_field, distribution_name, data_found)
            call check_data_found(data_field, data_found)
            select case (distribution_name)
                case ("rectangular")
                    allocate(Rectangular_Half_Distribution :: half_distribution)
                case default
                    call error_exit("procedures_half_distribution_factory: create: "//&
                        distribution_name//" distribution_name is unknown. Choose: 'rectangular'.")
            end select
        else
            allocate(Null_Half_Distribution :: half_distribution)
        end if

        select type (half_distribution)
            type is (Rectangular_Half_Distribution)
                data_field = prefix//"Half Distribution.delta distance"
                call exploring_data%get(data_field, delta_distance, data_found)
                call check_data_found(data_field, data_found)
                call half_distribution%set(delta_distance)
            type is (Null_Half_Distribution)
            class default
                call error_exit("procedures_half_distribution_factory: create: "//&
                        "distribution_name: unknown type.")
        end select
    end subroutine create

    subroutine destroy(half_distribution)
        class(Abstract_Half_Distribution), allocatable, intent(inout) :: half_distribution

        if (allocated(half_distribution)) deallocate(half_distribution)
    end subroutine destroy

end module procedures_half_distribution_factory
