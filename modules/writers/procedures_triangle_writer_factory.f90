module procedures_triangle_writer_factory

use types_component_wrapper, only: Component_Wrapper
use types_pair_potential_wrapper, only: Pair_Potentials_Line
use types_des_real_pair_wrapper, only: DES_Real_Pairs_Line
use classes_triangle_writer, only: Selectors_Line, &
    Abstract_Triangle_Writer, Concrete_Triangle_Writer, Null_Triangle_Writer
use procedures_property_inquirers, only: component_exists, components_interact

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_switches
    module procedure :: create_short_energies
    module procedure :: create_dipolar_energies
end interface create

contains

    subroutine create_switches(switches, components, filename)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: switches
        type(Component_Wrapper), intent(in) :: components(:)
        character(len=*), intent(in) :: filename

        type(Selectors_Line) :: selectors(size(components))
        logical :: some_components_exist, exist_ij
        integer :: i_component, j_component

        some_components_exist = .false.
        do j_component = 1, size(selectors)
            allocate(selectors(j_component)%line(j_component))
            do i_component = 1, size(selectors(j_component)%line)
                exist_ij = component_exists(components(i_component)%number) .and. &
                    component_exists(components(j_component)%number)
                some_components_exist = some_components_exist .or. exist_ij
                selectors(j_component)%line(i_component) = exist_ij .and. i_component /= j_component
            end do
        end do
        if (some_components_exist) then
            allocate(Concrete_Triangle_Writer :: switches)
        else
            allocate(Null_Triangle_Writer :: switches)
        end if
        call switches%construct(selectors, filename)
        call deallocate_selectors(selectors)
    end subroutine create_switches

    subroutine create_short_energies(energies, pairs, filename)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: energies
        type(Pair_Potentials_Line), intent(in) :: pairs(:)
        character(len=*), intent(in) :: filename

        type(Selectors_Line) :: selectors(size(pairs))
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
        call deallocate_selectors(selectors)
    end subroutine create_short_energies

    subroutine create_dipolar_energies(energies, pairs, filename)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: energies
        type(DES_Real_Pairs_Line), intent(in) :: pairs(:)
        character(len=*), intent(in) :: filename

        type(Selectors_Line) :: selectors(size(pairs))
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
        call deallocate_selectors(selectors)
    end subroutine create_dipolar_energies

    subroutine deallocate_selectors(selectors)
        type(Selectors_Line), intent(inout) :: selectors(:)

        integer :: i_component
        do i_component = size(selectors), 1, -1
            if (allocated(selectors(i_component)%line)) then
                deallocate(selectors(i_component)%line)
            end if
        end do
    end subroutine deallocate_selectors

    subroutine destroy(triangle)
        class(Abstract_Triangle_Writer), allocatable, intent(inout) :: triangle

        if (allocated(triangle)) then
            call triangle%destroy()
            deallocate(triangle)
        end if
    end subroutine destroy

end module procedures_triangle_writer_factory
