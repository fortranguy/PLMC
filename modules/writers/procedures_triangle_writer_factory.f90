module procedures_triangle_writer_factory

use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only: component_exists, component_can_exchange
use types_pair_potential_wrapper, only: Pair_Potentials_Line
use procedures_short_interactions_inquirers, only: components_interact
use types_logical_line, only: Concrete_Logical_Line
use classes_triangle_writer, only: Abstract_Triangle_Writer, Concrete_Triangle_Writer, &
    Null_Triangle_Writer

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_short_energies
    module procedure :: create_dipolar_energies
    module procedure :: create_switches
end interface create

contains

    subroutine create_short_energies(energies, pairs, filename)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: energies
        type(Pair_Potentials_Line), intent(in) :: pairs(:)
        character(len=*), intent(in) :: filename

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

        if (some_components_interact) then
            allocate(Concrete_Triangle_Writer :: energies)
        else
            allocate(Null_Triangle_Writer :: energies)
        end if
        call energies%construct(selectors, filename)
    end subroutine create_short_energies

    subroutine create_dipolar_energies(energies, are_dipolar, filename)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: energies
        logical, intent(in) :: are_dipolar(:)
        character(len=*), intent(in) :: filename

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

        if (some_components_interact) then
            allocate(Concrete_Triangle_Writer :: energies)
        else
            allocate(Null_Triangle_Writer :: energies)
        end if
        call energies%construct(selectors, filename)
    end subroutine create_dipolar_energies

    subroutine create_switches(switches, components, filename)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: switches
        type(Component_Wrapper), intent(in) :: components(:)
        character(len=*), intent(in) :: filename

        type(Concrete_Logical_Line) :: selectors(size(components))
        logical :: some_couples_exist, exist_ij
        integer :: i_component, j_component

        some_couples_exist = .false.
        do j_component = 1, size(selectors)
            allocate(selectors(j_component)%line(j_component))
            do i_component = 1, size(selectors(j_component)%line)
                exist_ij = &
                    component_exists(components(i_component)%number) .and. &
                    component_exists(components(j_component)%number) .and. &
                    i_component /= j_component
                some_couples_exist = some_couples_exist .or. exist_ij
                selectors(j_component)%line(i_component) = exist_ij
            end do
        end do
        if (some_couples_exist) then
            allocate(Concrete_Triangle_Writer :: switches)
        else
            allocate(Null_Triangle_Writer :: switches)
        end if
        call switches%construct(selectors, filename)
    end subroutine create_switches

    subroutine destroy(triangle)
        class(Abstract_Triangle_Writer), allocatable, intent(inout) :: triangle

        if (allocated(triangle)) then
            call triangle%destroy()
            deallocate(triangle)
        end if
    end subroutine destroy

end module procedures_triangle_writer_factory
