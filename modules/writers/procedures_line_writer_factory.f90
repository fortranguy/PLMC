module procedures_line_writer_factory

use classes_external_field, only: Abstract_External_Field
use types_component_wrapper, only: Component_Wrapper
use types_pair_potential_wrapper, only: Pair_Potential_Wrapper
use classes_line_writer, only: Abstract_Line_Writer, Concrete_Line_Writer, Null_Line_Writer
use procedures_property_inquirers, only: apply_external_field, component_interacts_with_wall, &
    component_is_dipolar

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_field
    module procedure :: create_walls
    module procedure :: create_line
end interface

contains

    subroutine create_field(field, external_field, components, filename)
        class(Abstract_Line_Writer), allocatable, intent(out) :: field
        class(Abstract_External_Field), intent(in) :: external_field
        type(Component_Wrapper), intent(in) :: components(:)
        character(len=*), intent(in) :: filename

        logical :: selector(size(components))
        integer :: i_component

        do i_component = 1, size(selector)
            selector(i_component) = component_is_dipolar(components(i_component)%dipolar_moments)
        end do

        if (any(selector) .and. apply_external_field(external_field)) then
            allocate(Concrete_Line_Writer :: field)
        else
            allocate(Null_Line_Writer :: field)
        end if
        call field%construct(selector, filename)
    end subroutine create_field

    subroutine create_walls(walls, wall_pairs, filename)
        class(Abstract_Line_Writer), allocatable, intent(out) :: walls
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        character(len=*), intent(in) :: filename

        logical :: selector(size(wall_pairs))
        integer :: i_component

        do i_component = 1, size(wall_pairs)
            selector(i_component) = component_interacts_with_wall(wall_pairs(i_component)%potential)
        end do

        if (any(selector)) then
            allocate(Concrete_Line_Writer :: walls)
        else
            allocate(Null_Line_Writer :: walls)
        end if
        call walls%construct(selector, filename)
    end subroutine create_walls

    subroutine create_line(line, selector, filename)
        class(Abstract_Line_Writer), allocatable, intent(out) :: line
        logical, intent(in) :: selector(:)
        character(len=*), intent(in) :: filename

        if (any(selector)) then
            allocate(Concrete_Line_Writer :: line)
        else
            allocate(Null_Line_Writer :: line)
        end if
        call line%construct(selector, filename)
    end subroutine create_line

    subroutine destroy(line)
        class(Abstract_Line_Writer), allocatable, intent(inout) :: line

        if (allocated(line)) then
            call line%destroy()
            deallocate(line)
        end if
    end subroutine destroy

end module procedures_line_writer_factory

