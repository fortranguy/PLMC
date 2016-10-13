module procedures_triangle_writer_factory

use types_logical_line, only: Concrete_Logical_Line
use types_string_wrapper, only: String_Wrapper
use procedures_environment_inquirers, only: box_size_can_change
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only: component_exists, component_can_exchange
use types_pair_potential_wrapper, only: Pair_Potentials_Line
use procedures_short_interactions_inquirers, only: components_interact
use types_changed_box_size_wrapper, only: Changed_Box_Size_Line
use classes_triangle_writer, only: Abstract_Triangle_Writer, Concrete_Triangle_Writer, &
    Null_Triangle_Writer

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_volumes_change_success
    module procedure :: create_short_energies
    module procedure :: create_dipolar_energies
    module procedure :: create_switches_successes
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_line
end interface destroy

contains

    subroutine create_volumes_change_success(successes, filename, changed_boxes_size)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: successes
        character(len=*), intent(in) :: filename
        type(Changed_Box_Size_Line), intent(in) :: changed_boxes_size(:)

        type(Concrete_Logical_Line) :: selectors(size(changed_boxes_size))
        logical :: some_boxes_can_change, can_change_ij
        integer :: i_box, j_box

        some_boxes_can_change = .false.
        do j_box = 1, size(selectors)
            allocate(selectors(j_box)%line(j_box))
            do i_box = 1, size(selectors(j_box)%line)
                can_change_ij = box_size_can_change(changed_boxes_size(j_box)%line(i_box)%changed)
                some_boxes_can_change = some_boxes_can_change .or. can_change_ij
                selectors(j_box)%line(i_box) = can_change_ij
            end do
        end do

        if (some_boxes_can_change) then
            allocate(Concrete_Triangle_Writer :: successes)
        else
            allocate(Null_Triangle_Writer :: successes)
        end if
        call successes%construct(filename, selectors)
    end subroutine create_volumes_change_success

    subroutine create_short_energies(energies, filename, pairs, visit_energies)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: energies
        character(len=*), intent(in) :: filename
        type(Pair_Potentials_Line), intent(in) :: pairs(:)
        logical, intent(in) :: visit_energies

        type(Concrete_Logical_Line) :: selectors(size(pairs))
        logical :: some_components_interact, interact_ij
        integer :: i_component, j_component

        some_components_interact = .false.
        do j_component = 1, size(selectors)
            allocate(selectors(j_component)%line(j_component))
            do i_component = 1, size(selectors(j_component)%line)
                interact_ij = components_interact(pairs(j_component)%line(i_component)%potential)
                some_components_interact = some_components_interact .or. interact_ij
                selectors(j_component)%line(i_component) = interact_ij
            end do
        end do

        if (some_components_interact .and. visit_energies) then
            allocate(Concrete_Triangle_Writer :: energies)
        else
            allocate(Null_Triangle_Writer :: energies)
        end if
        call energies%construct(filename, selectors)
    end subroutine create_short_energies

    subroutine create_dipolar_energies(energies, filename, are_dipolar, visit_energies)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: energies
        character(len=*), intent(in) :: filename
        logical, intent(in) :: are_dipolar(:)
        logical, intent(in) :: visit_energies

        type(Concrete_Logical_Line) :: selectors(size(are_dipolar))
        logical :: some_components_interact, interact_ij
        integer :: i_component, j_component

        some_components_interact = .false.
        do j_component = 1, size(selectors)
            allocate(selectors(j_component)%line(j_component))
            do i_component = 1, size(selectors(j_component)%line)
                interact_ij = are_dipolar(i_component) .and. are_dipolar(j_component)
                some_components_interact = some_components_interact .or. interact_ij
                selectors(j_component)%line(i_component) = interact_ij
            end do
        end do

        if (some_components_interact .and. visit_energies) then
            allocate(Concrete_Triangle_Writer :: energies)
        else
            allocate(Null_Triangle_Writer :: energies)
        end if
        call energies%construct(filename, selectors)
    end subroutine create_dipolar_energies

    subroutine destroy_element(writer)
        class(Abstract_Triangle_Writer), allocatable, intent(inout) :: writer

        if (allocated(writer)) then
            call writer%destroy()
            deallocate(writer)
        end if
    end subroutine destroy_element

    subroutine create_switches_successes(successes, paths, filename, components)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: successes(:)
        type(String_Wrapper), intent(in) :: paths(:)
        character(len=*), intent(in) :: filename
        type(Component_Wrapper), intent(in) :: components(:, :)

        type(Concrete_Logical_Line) :: selectors(size(components, 1), size(components, 2))
        logical :: some_couples_exist(size(selectors, 2)), exist_ij
        integer :: i_box, i_component, j_component

        some_couples_exist = .false.
        do i_box = 1, size(selectors, 2)
            do j_component = 1, size(selectors, 1)
                allocate(selectors(j_component, i_box)%line(j_component))
                do i_component = 1, size(selectors(j_component, i_box)%line)
                    exist_ij = component_exists(components(i_component, i_box)%num_particles) .and.&
                        component_exists(components(j_component, i_box)%num_particles) .and. &
                        i_component /= j_component
                    some_couples_exist(i_box) = some_couples_exist(i_box) .or. exist_ij
                    selectors(j_component, i_box)%line(i_component) = exist_ij
                end do
            end do
        end do

        if (all(some_couples_exist)) then
            allocate(Concrete_Triangle_Writer :: successes(size(paths)))
        else
            allocate(Null_Triangle_Writer :: successes(size(paths)))
        end if

        do i_box = 1, size(successes)
            call successes(i_box)%construct(paths(i_box)%string//filename, selectors(:, i_box))
        end do
    end subroutine create_switches_successes

    subroutine destroy_line(writers)
        class(Abstract_Triangle_Writer), allocatable, intent(inout) :: writers(:)

        integer :: i_element

        if (allocated(writers)) then
            do i_element = size(writers), 1, -1
                call writers(i_element)%destroy()
            end do
            deallocate(writers)
        end if
    end subroutine destroy_line

end module procedures_triangle_writer_factory
