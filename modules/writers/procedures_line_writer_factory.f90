module procedures_line_writer_factory

use classes_external_field, only: Abstract_External_Field
use procedures_environment_inquirers, only: apply_external_field
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only: component_is_dipolar
use types_pair_potential_wrapper, only: Pair_Potential_Wrapper
use procedures_short_interactions_inquirers, only: component_interacts_with_wall
use classes_line_writer, only: Abstract_Line_Writer, Concrete_Line_Writer, Null_Line_Writer

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_walls
    module procedure :: create_field
    module procedure :: create_line
end interface

contains

    subroutine create_walls(walls_energies, wall_pairs, visit_energies, filename)
        class(Abstract_Line_Writer), allocatable, intent(out) :: walls_energies
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        logical, intent(in) :: visit_energies
        character(len=*), intent(in) :: filename

        logical :: selector(size(wall_pairs))
        integer :: i_component

        do i_component = 1, size(wall_pairs)
            selector(i_component) = component_interacts_with_wall(wall_pairs(i_component)%potential)
        end do

        if (any(selector) .and. visit_energies) then
            allocate(Concrete_Line_Writer :: walls_energies)
        else
            allocate(Null_Line_Writer :: walls_energies)
        end if
        call walls_energies%construct(selector, filename)
    end subroutine create_walls

    subroutine create_field(field_energies, external_field, components, visit_energies, filename)
        class(Abstract_Line_Writer), allocatable, intent(out) :: field_energies
        class(Abstract_External_Field), intent(in) :: external_field
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: visit_energies
        character(len=*), intent(in) :: filename

        logical :: selector(size(components))
        integer :: i_component

        do i_component = 1, size(selector)
            selector(i_component) = component_is_dipolar(components(i_component)%dipole_moments)
        end do

        if (any(selector) .and. apply_external_field(external_field) .and. visit_energies) then
            allocate(Concrete_Line_Writer :: field_energies)
        else
            allocate(Null_Line_Writer :: field_energies)
        end if
        call field_energies%construct(selector, filename)
    end subroutine create_field

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

