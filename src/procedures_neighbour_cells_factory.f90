module procedures_neighbour_cells_factory

use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_environment_inquirers, only: periodicity_is_xyz, periodicity_is_xy
use classes_pair_potential, only: Abstract_Pair_Potential, Pair_Potential_Line
use classes_hard_contact, only: Abstract_Hard_Contact
use classes_dipoles_neighbourhood ,only: Abstract_Dipolar_Neighbourhood, Dipolar_Neighbourhood_Line
use classes_neighbour_cells, only: Abstract_Neighbour_Cells, XYZ_PBC_Neighbour_Cells, &
    XY_PBC_Neighbour_Cells, Null_Neighbour_Cells, Neighbour_Cells_Line

implicit none

private
public :: create_triangle, create_element, destroy_triangle, destroy_element, allocate_triangle

contains

    subroutine create_triangle(cells, periodic_box, accessible_domain, components_pairs, &
        hard_contact, dipolar_neighbourhoods, interact)
        type(Neighbour_Cells_Line), allocatable, intent(out) :: cells(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Parallelepiped_Domain), intent(in) :: accessible_domain
        type(Pair_Potential_Line), intent(in) :: components_pairs(:)
        class(Abstract_Hard_Contact), intent(in) :: hard_contact
        type(Dipolar_Neighbourhood_Line), intent(in) :: dipolar_neighbourhoods(:)
        logical, intent(in) :: interact

        integer :: i_component, j_component

        call allocate_triangle(cells, size(components_pairs))
        do j_component = 1, size(cells)
            do i_component = 1, size(cells(j_component)%line)
                associate(pair_ij => components_pairs(j_component)%line(i_component)%potential, &
                        neighbourhood_ij => dipolar_neighbourhoods(j_component)%line(i_component)%&
                        neighbourhood)
                    call create_element(cells(j_component)%line(i_component)%cells, periodic_box, &
                        accessible_domain, pair_ij, hard_contact, neighbourhood_ij, interact)
                end associate
            end do
        end do
    end subroutine create_triangle

    subroutine create_element(cells, periodic_box, accessible_domain, pair_potential, hard_contact,&
        dipolar_neighbourhood, interact)
        class(Abstract_Neighbour_Cells), allocatable, intent(out) :: cells
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Parallelepiped_Domain), intent(in) :: accessible_domain
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        class(Abstract_Hard_Contact), intent(in) :: hard_contact
        class(Abstract_Dipolar_Neighbourhood), intent(in) :: dipolar_neighbourhood
        logical, intent(in) :: interact

        if (interact) then
            if (periodicity_is_xyz(periodic_box)) then
                allocate(XYZ_PBC_Neighbour_Cells :: cells)
            else if (periodicity_is_xy(periodic_box)) then
                allocate(XY_PBC_Neighbour_Cells :: cells)
            else
                call error_exit("procedures_neighbour_cells_factory: create_element: "//&
                        "box periodicity is unknown.")
            end if
        else
            allocate(Null_Neighbour_Cells :: cells)
        end if
        call cells%construct(accessible_domain, pair_potential, hard_contact, dipolar_neighbourhood)
    end subroutine create_element

    subroutine allocate_triangle(cells, num_components)
        type(Neighbour_Cells_Line), allocatable, intent(out) :: cells(:)
        integer, intent(in) :: num_components

        integer :: i_component

        allocate(cells(num_components))
        do i_component = 1, size(cells)
            allocate(cells(i_component)%line(i_component))
        end do
    end subroutine allocate_triangle

    subroutine destroy_triangle(cells)
        type(Neighbour_Cells_Line), allocatable, intent(inout) :: cells(:)

        integer :: i_component, j_component

        if (allocated(cells)) then
            do j_component = size(cells), 1, -1
                if (allocated(cells(j_component)%line)) then
                    do i_component = size(cells(j_component)%line), 1, -1
                        call destroy_element(cells(j_component)%line(i_component)%cells)
                    end do
                    deallocate(cells(j_component)%line)
                end if
            end do
            deallocate(cells)
        end if
    end subroutine destroy_triangle

    subroutine destroy_element(cells)
        class(Abstract_Neighbour_Cells), allocatable, intent(inout) :: cells

        if (allocated(cells)) then
            call cells%destroy()
            deallocate(cells)
        end if
    end subroutine destroy_element

end module procedures_neighbour_cells_factory
