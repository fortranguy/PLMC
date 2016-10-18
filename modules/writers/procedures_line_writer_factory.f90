module procedures_line_writer_factory

use types_string_wrapper, only: String_Wrapper
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
    module procedure :: create_walls_energies
    module procedure :: create_field_energies
    module procedure :: create_line
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_line
end interface

contains

    subroutine create_walls_energies(energies, filename, wall_pairs, visit_energies)
        class(Abstract_Line_Writer), allocatable, intent(out) :: energies
        character(len=*), intent(in) :: filename
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        logical, intent(in) :: visit_energies

        logical :: selector(size(wall_pairs))
        integer :: i_component

        do i_component = 1, size(wall_pairs)
            selector(i_component) = component_interacts_with_wall(wall_pairs(i_component)%potential)
        end do

        if (any(selector) .and. visit_energies) then
            allocate(Concrete_Line_Writer :: energies)
        else
            allocate(Null_Line_Writer :: energies)
        end if
        call energies%construct(filename, selector)
    end subroutine create_walls_energies

    subroutine create_field_energies(energies, filename, external_field, components, visit_energies)
        class(Abstract_Line_Writer), allocatable, intent(out) :: energies
        character(len=*), intent(in) :: filename
        class(Abstract_External_Field), intent(in) :: external_field
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: visit_energies

        logical :: selector(size(components))
        integer :: i_component

        do i_component = 1, size(selector)
            selector(i_component) = component_is_dipolar(components(i_component)%dipole_moments)
        end do

        if (any(selector) .and. apply_external_field(external_field) .and. visit_energies) then
            allocate(Concrete_Line_Writer :: energies)
        else
            allocate(Null_Line_Writer :: energies)
        end if
        call energies%construct(filename, selector)
    end subroutine create_field_energies

    !> @note Condition may not be optimal.
    subroutine create_line(writers, paths, filename, selectors)
        class(Abstract_Line_Writer), allocatable, intent(out) :: writers(:)
        type(String_Wrapper), intent(in) :: paths(:)
        character(len=*), intent(in) :: filename
        logical, intent(in) :: selectors(:, :)

        integer :: i_element

        if (any(selectors)) then
            allocate(Concrete_Line_Writer :: writers(size(paths)))
        else
            allocate(Null_Line_Writer :: writers(size(paths)))
        end if

        do i_element = 1, size(writers)
            call writers(i_element)%construct(paths(i_element)%string//filename, &
                selectors(:, i_element))
        end do
    end subroutine create_line

    subroutine destroy_line(writers)
        class(Abstract_Line_Writer), allocatable, intent(inout) :: writers(:)

        integer :: i_element

        if (allocated(writers)) then
            do i_element = size(writers), 1, -1
                call writers(i_element)%destroy()
            end do
            deallocate(writers)
        end if
    end subroutine destroy_line

    subroutine destroy_element(writer)
        class(Abstract_Line_Writer), allocatable, intent(inout) :: writer

        if (allocated(writer)) then
            call writer%destroy()
            deallocate(writer)
        end if
    end subroutine destroy_element

end module procedures_line_writer_factory

