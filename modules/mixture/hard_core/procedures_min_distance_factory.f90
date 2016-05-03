module procedures_min_distance_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_min_distance, only: Abstract_Min_Distance, Concrete_Min_Distance, Null_Min_Distance

implicit none

private
public :: create, destroy

contains

    subroutine create(min_distance, exists, generating_data, prefix)
        class(Abstract_Min_Distance), allocatable, intent(out) :: min_distance
        logical, intent(in) :: exists
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: min_distance_value

        if (exists) then
            allocate(Concrete_Min_Distance :: min_distance)
            data_field = prefix//"minimum distance"
            call generating_data%get(data_field, min_distance_value, data_found)
            call check_data_found(data_field, data_found)
        else
            allocate(Null_Min_Distance :: min_distance)
        end if
        call min_distance%set(min_distance_value)
    end subroutine create

    subroutine destroy(min_distance)
        class(Abstract_Min_Distance), allocatable, intent(inout) :: min_distance

        if (allocated(min_distance)) deallocate(min_distance)
    end subroutine destroy

end module procedures_min_distance_factory
