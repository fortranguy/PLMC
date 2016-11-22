module procedures_triangle_writer_factory

use types_logical_wrapper, only: Logical_Line
use types_string_wrapper, only: String_Wrapper
use procedures_environment_inquirers, only: box_size_can_change
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_properties, only: set_have_positions
use classes_pair_potential, only: Pair_Potential_Line
use procedures_short_interactions_inquirers, only: components_interact
use classes_exchanged_boxes_size, only: Exchanged_Boxes_Size_Line
use classes_triangle_writer, only: Abstract_Triangle_Writer, Concrete_Triangle_Writer, &
    Null_Triangle_Writer

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_volumes_exchange_success
    module procedure :: create_short_energies
    module procedure :: create_dipolar_energies
    module procedure :: create_switches_successes
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_line
end interface destroy

contains

    subroutine create_volumes_exchange_success(successes, filename, exchanged_boxes_size)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: successes
        character(len=*), intent(in) :: filename
        type(Exchanged_Boxes_Size_Line), intent(in) ::exchanged_boxes_size(:)

        type(Logical_Line) :: selectors(size(exchanged_boxes_size))
        logical :: some_boxes_exchange_volume, exchange_ij
        integer :: i_box, j_box

        some_boxes_exchange_volume = .false.
        do j_box = 1, size(selectors)
            allocate(selectors(j_box)%line(j_box))
            do i_box = 1, size(selectors(j_box)%line)
                exchange_ij = i_box /= j_box
                some_boxes_exchange_volume = some_boxes_exchange_volume .or. exchange_ij
                selectors(j_box)%line(i_box) = exchange_ij
            end do
        end do

        if (some_boxes_exchange_volume) then
            allocate(Concrete_Triangle_Writer :: successes)
        else
            allocate(Null_Triangle_Writer :: successes)
        end if
        call successes%construct(filename, selectors)
    end subroutine create_volumes_exchange_success

    subroutine create_short_energies(energies, filename, pairs, visit_energies)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: energies
        character(len=*), intent(in) :: filename
        type(Pair_Potential_Line), intent(in) :: pairs(:)
        logical, intent(in) :: visit_energies

        type(Logical_Line) :: selectors(size(pairs))
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

        type(Logical_Line) :: selectors(size(are_dipolar))
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

        type(Logical_Line) :: selectors(size(components, 1), size(components, 2))
        logical :: have_positions(size(components, 1), size(components, 2)), have_positions_ij
        logical :: some_couples_have_positions(size(selectors, 2))
        integer :: i_box, i_component, j_component

        call set_have_positions(have_positions, components)
        some_couples_have_positions = .false.
        do i_box = 1, size(selectors, 2)
            do j_component = 1, size(selectors, 1)
                allocate(selectors(j_component, i_box)%line(j_component))
                do i_component = 1, size(selectors(j_component, i_box)%line)
                    have_positions_ij = have_positions(i_component, i_box) .and. &
                        have_positions(j_component, i_box) .and. i_component /= j_component
                    some_couples_have_positions(i_box) = some_couples_have_positions(i_box) .or. &
                        have_positions_ij
                    selectors(j_component, i_box)%line(i_component) = have_positions_ij
                end do
            end do
        end do

        if (all(some_couples_have_positions)) then
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
