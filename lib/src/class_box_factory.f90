module class_box_factory

use json_module, only: json_file
use module_data, only: test_data_found
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use types_box, only: Box_Wrapper

implicit none

private

    type, public :: Concrete_Box_Factory
    private
        type(json_file), pointer :: input_data
        character(len=:), allocatable :: prefix
    contains
        procedure :: allocate => Concrete_Box_Factory_allocate
        procedure :: allocate_periodic_box => Concrete_Box_Factory_allocate_periodic_box
        !procedure :: construct => Concrete_Box_Factory_construct
        procedure :: destroy => Concrete_Box_Factory_destroy
    end type Concrete_Box_Factory

contains

    subroutine Concrete_Box_Factory_allocate(this, box, input_data, prefix)
        class(Concrete_Box_Factory), intent(out) :: this
        type(Box_Wrapper), intent(out) :: box
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        this%input_data = input_data
        this%prefix = prefix
        call this%allocate_periodic_box(box%periodic_box)
    end subroutine Concrete_Box_Factory_allocate

    subroutine Concrete_Box_Factory_allocate_periodic_box(this, periodic_box)
        class(Concrete_Box_Factory), intent(in) :: this
        class(Abstract_Periodic_Box), allocatable, intent(out) :: periodic_box

        character(len=:), allocatable :: data_field
        logical :: data_found
        character(len=:), allocatable :: box_name

        data_field = "Periodic Box.name"
        call this%input_data%get(data_field, box_name, data_found)
        call test_data_found(data_field, data_found)
        select case(box_name)
            case("XYZ")
                allocate(XYZ_Periodic_Box :: periodic_box)
            case("XY")
                allocate(XY_Periodic_Box :: periodic_box)
            case default
                call error_exit(data_field//" unkown.")
        end select
        deallocate(box_name)
    end subroutine Concrete_Box_Factory_allocate_periodic_box

    subroutine Concrete_Box_Factory_destroy(this, box)
        class(Concrete_Box_Factory), intent(inout) :: this
        type(Box_Wrapper), intent(inout) :: box

        if (allocated(box%periodic_box)) deallocate(box%periodic_box)
        if (allocated(this%prefix)) deallocate(this%prefix)
        this%input_data => null()
    end subroutine Concrete_Box_Factory_destroy

end module class_box_factory
