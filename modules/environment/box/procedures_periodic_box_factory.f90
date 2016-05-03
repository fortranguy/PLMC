module procedures_periodic_box_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box

implicit none

private
public :: create, destroy

contains

    subroutine create(periodic_box, generating_data, prefix)
        class(Abstract_Periodic_Box), allocatable, intent(out) :: periodic_box
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: box_periodicity
        real(DP), allocatable :: box_size(:)
        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = prefix//"Box.periodicity"
        call generating_data%get(data_field, box_periodicity, data_found)
        call check_data_found(data_field, data_found)
        select case (box_periodicity)
            case ("XYZ")
                allocate(XYZ_Periodic_Box :: periodic_box)
            case ("XY")
                allocate(XY_Periodic_Box :: periodic_box)
            case default
                call error_exit(data_field//" unknown. Choose between: 'XYZ' and 'XY'")
        end select
        data_field = prefix//"Box.initial size"
        call generating_data%get(data_field, box_size, data_found)
        call check_data_found(data_field, data_found)
        call periodic_box%set(box_size)
    end subroutine create

    subroutine destroy(periodic_box)
        class(Abstract_Periodic_Box), allocatable, intent(inout) :: periodic_box

        if (allocated(periodic_box)) deallocate(periodic_box)
    end subroutine destroy

end module procedures_periodic_box_factory
