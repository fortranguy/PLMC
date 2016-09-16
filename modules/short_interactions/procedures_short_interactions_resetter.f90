module procedures_short_interactions_resetter

use types_logical_line, only: Concrete_Logical_Line
use types_neighbour_cells_wrapper, only: Neighbour_Cells_Line
use classes_visitable_cells, only: Abstract_Visitable_Cells

implicit none

private
public :: reset

interface reset
    module procedure :: reset_neighbour_cells
    module procedure :: reset_visitable_cells
    module procedure :: box_size_change_reset_cells
end interface reset

contains

    subroutine reset_neighbour_cells(cells)
        type(Neighbour_Cells_Line), intent(inout) :: cells(:)

        integer :: i_component, j_component
        logical :: only_resized_dummy

        do j_component = 1, size(cells)
            do i_component = 1, size(cells(j_component)%line)
                call cells(j_component)%line(i_component)%cells%reset(only_resized_dummy)
            end do
        end do
    end subroutine reset_neighbour_cells

    subroutine reset_visitable_cells(cells)
        class(Abstract_Visitable_Cells), intent(inout) :: cells(:, :)

        integer :: i_component, j_component

        do j_component = 1, size(cells, 2)
            do i_component = 1, size(cells, 1)
                call cells(i_component, j_component)%reset()
            end do
        end do
    end subroutine reset_visitable_cells

    subroutine box_size_change_reset_cells(neighbour_cells, only_resized_triangle, visitable_cells)
        type(Neighbour_Cells_Line), intent(inout) :: neighbour_cells(:)
        type(Concrete_Logical_Line), allocatable, intent(out) :: only_resized_triangle(:)
        class(Abstract_Visitable_Cells), intent(inout) :: visitable_cells(:, :)

        integer :: i_component, j_component

        allocate(only_resized_triangle(size(neighbour_cells)))
        do j_component = 1, size(neighbour_cells)
            allocate(only_resized_triangle(j_component)%line(j_component)) ! allocate elsewhere?
            do i_component = 1, size(neighbour_cells(j_component)%line)
                call neighbour_cells(j_component)%line(i_component)%cells%&
                    reset(only_resized_triangle(j_component)%line(i_component))
                if (.not. only_resized_triangle(j_component)%line(i_component)) then
                    call visitable_cells(i_component, j_component)%reset()
                    if (i_component /= j_component) call visitable_cells(j_component, i_component)%&
                        reset()
                end if
            end do
        end do
    end subroutine box_size_change_reset_cells

end module procedures_short_interactions_resetter
