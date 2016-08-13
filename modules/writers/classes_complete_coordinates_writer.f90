module classes_complete_coordinates_writer

use classes_box_size_writer, only: Abstract_Box_Size_Writer
use types_component_coordinates_writer_wrapper, only: Component_Coordinates_Writer_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Complete_Coordinates_Writer
    private
        class(Abstract_Box_Size_Writer), allocatable :: box_size
        type(Component_Coordinates_Writer_Wrapper), allocatable :: components_coordinates(:)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
    end type Abstract_Complete_Coordinates_Writer

contains

!implementation Abstract_Complete_Coordinates_Writer

    subroutine Abstract_construct(this, box_size, components_coordinates)
        class(Abstract_Complete_Coordinates_Writer), intent(out) :: this
        class(Abstract_Box_Size_Writer), intent(in) :: box_size
        type(Component_Coordinates_Writer_Wrapper), intent(in) :: components_coordinates(:)

        allocate(this%box_size, source=box_size)
        allocate(this%components_coordinates, source=components_coordinates)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Complete_Coordinates_Writer), intent(inout) :: this

        integer :: i_component

        if (allocated(this%box_size)) then
            call this%box_size%destroy()
            deallocate(this%box_size)
        end if
        if (allocated(this%components_coordinates)) then !factory?
            do i_component = size(this%components_coordinates), 1, -1
                if (allocated(this%components_coordinates(i_component)%writer)) then
                    call this%components_coordinates(i_component)%writer%destroy()
                    deallocate(this%components_coordinates(i_component)%writer)
                end if
            end do
        end if
    end subroutine Abstract_destroy

!end implementation Abstract_Complete_Coordinates_Writer

end module classes_complete_coordinates_writer
