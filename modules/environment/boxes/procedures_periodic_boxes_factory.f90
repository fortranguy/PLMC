module procedures_periodic_boxes_factory

use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use procedures_environment_inquirers, only: property_num_boxes => num_boxes

implicit none

private
public :: create, destroy

contains

    subroutine create(periodic_boxes, generating_data, prefix)
        class(Abstract_Periodic_Box), allocatable, intent(out) :: periodic_boxes(:)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: num_boxes
        character(len=:), allocatable :: boxes_periodicity
        character(len=:), allocatable :: data_field
        logical :: data_found

        num_boxes = property_num_boxes(generating_data, prefix)
        data_field = prefix//"Boxes.periodicity"
        call generating_data%get(data_field, boxes_periodicity, data_found)
        call check_data_found(data_field, data_found)
        select case (boxes_periodicity)
            case ("XYZ")
                allocate(XYZ_Periodic_Box :: periodic_boxes(num_boxes))
            case ("XY")
                allocate(XY_Periodic_Box :: periodic_boxes(num_boxes))
            case default
                call error_exit(data_field//" unknown. Choose between: 'XYZ' and 'XY'")
        end select
    end subroutine create

    subroutine destroy(periodic_boxes)
        class(Abstract_Periodic_Box), allocatable, intent(inout) :: periodic_boxes(:)

        if (allocated(periodic_boxes)) deallocate(periodic_boxes)
    end subroutine destroy

end module procedures_periodic_boxes_factory
