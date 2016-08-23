module procedures_visitable_cells_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use types_component_wrapper, only: Component_Wrapper
use classes_hard_contact, only: Abstract_Hard_Contact
use types_pair_potential_wrapper, only: Pair_Potentials_Line
use types_neighbour_cells_wrapper, only: Neighbour_Cells_Line
use classes_visitable_list, only: Abstract_Visitable_List
use classes_visitable_cells, only: Abstract_Visitable_Cells, Concrete_Visitable_Cells, &
    Null_Visitable_Cells

implicit none

private
public :: create, destroy

contains

    subroutine create(cells, periodic_box, components, hard_contact, pairs, neighbour_cells, &
        list_mold, interact)
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: cells(:, :)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: components(:)
        class(Abstract_Hard_Contact), intent(in) :: hard_contact
        type(Pair_Potentials_Line), intent(in) :: pairs(:)
        type(Neighbour_Cells_Line), intent(in) :: neighbour_cells(:)
        class(Abstract_Visitable_List), intent(in) :: list_mold
        logical, intent(in) :: interact

        integer :: i_component, j_component
        integer :: i_pair, j_pair

        if (interact) then
            allocate(Concrete_Visitable_Cells :: cells(size(pairs), size(pairs)))
        else
            allocate(Null_Visitable_Cells :: cells(size(pairs), size(pairs)))
        end if

        do j_component = 1, size(cells, 2)
            do i_component = 1, size(cells, 1)
                j_pair = maxval([i_component, j_component])
                i_pair = minval([i_component, j_component])
                associate (pair_ij => pairs(j_pair)%line(i_pair)%potential, &
                    neighbours_ij => neighbour_cells(j_pair)%line(i_pair)%cells)
                    call cells(i_component, j_component)%construct(periodic_box, &
                        components(i_component)%positions, hard_contact, pair_ij, list_mold)
                    call cells(i_component, j_component)%set(neighbours_ij)
                end associate
            end do
        end do
    end subroutine create

    subroutine destroy(cells)
        class(Abstract_Visitable_Cells), allocatable, intent(inout) :: cells(:, :)

        integer :: i_component, j_component

        if (allocated(cells)) then
            do j_component = size(cells, 2), 1, -1
                do i_component = size(cells, 1), 1, -1
                    call cells(i_component, j_component)%destroy()
                end do
            end do
            deallocate(cells)
        end if
    end subroutine destroy

end module procedures_visitable_cells_factory
