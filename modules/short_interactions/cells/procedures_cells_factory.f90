module procedures_cells_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_hard_contact, only: Abstract_Hard_Contact
use types_pair_potential_wrapper, only: Pair_Potentials_Line
use types_neighbour_cells_wrapper, only: Neighbour_Cells_Line
use procedures_neighbour_cells_factory, only: neighbour_cells_create => create, &
    neighbour_cells_destroy => destroy
use procedures_visitable_cells_factory, only: visitable_cells_create => create, &
    visitable_cells_destroy => destroy

implicit none

private
public :: create, destroy

interface create
    module procedure :: visitable_cells_create
    module procedure :: create_components
    module procedure :: neighbour_cells_create
end interface create

interface destroy
    module procedure :: neighbour_cells_destroy
    module procedure :: destroy_components
    module procedure :: visitable_cells_destroy
end interface destroy

contains

    subroutine create_components(cells, periodic_box, accessible_domain, hard_contact, pairs, &
        interact)
        type(Neighbour_Cells_Line), allocatable, intent(out) :: cells(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Parallelepiped_Domain), intent(in) :: accessible_domain
        class(Abstract_Hard_Contact), intent(in) :: hard_contact
        type(Pair_Potentials_Line), intent(in) :: pairs(:)
        logical, intent(in) :: interact

        integer :: i_component, j_component

        allocate(cells(size(pairs)))
        do j_component = 1, size(cells)
            allocate(cells(j_component)%line(j_component))
            do i_component = 1, size(cells(j_component)%line)
                associate(pair_ij => pairs(j_component)%line(i_component)%potential)
                    call create(cells(j_component)%line(i_component)%cells, periodic_box, &
                        accessible_domain, hard_contact, pair_ij, interact)
                end associate
            end do
        end do
    end subroutine create_components

    subroutine destroy_components(cells)
        type(Neighbour_Cells_Line), allocatable, intent(inout) :: cells(:)

        integer :: i_component, j_component

        if (allocated(cells)) then
            do j_component = size(cells), 1, -1
                if (allocated(cells(j_component)%line)) then
                    do i_component = size(cells(j_component)%line), 1, -1
                        call destroy(cells(j_component)%line(i_component)%cells)
                    end do
                    deallocate(cells(j_component)%line)
                end if
            end do
            deallocate(cells)
        end if
    end subroutine destroy_components

end module procedures_cells_factory
