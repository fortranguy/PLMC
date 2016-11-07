module procedures_cells_memento

use types_logical_wrapper, only: Logical_Line
use classes_visitable_cells_memento, only: Abstract_Visitable_Cells_Memento
use types_cells_wrapper, only: Cells_Wrapper
use procedures_cells_factory, only: cells_destroy => destroy, cells_allocate_triangle => &
    allocate_triangle

implicit none

private
public :: save, restore

contains

    subroutine save(cells_target, only_resized_triangle, visitable_cells_memento, cells_source)
        type(Cells_Wrapper), intent(out) :: cells_target
        type(Logical_Line), intent(inout) ::only_resized_triangle(:)
        class(Abstract_Visitable_Cells_Memento), intent(in) :: visitable_cells_memento
        type(Cells_Wrapper), intent(in) :: cells_source

        integer :: i_component, j_component

        call cells_allocate_triangle(cells_target%neighbour_cells, &
            size(cells_source%neighbour_cells))
        do j_component = 1, size(cells_target%neighbour_cells)
            do i_component = 1, size(cells_target%neighbour_cells(j_component)%line)
                only_resized_triangle(j_component)%line(i_component) = cells_source%&
                    neighbour_cells(j_component)%line(i_component)%cells%resize_only()
                if (.not. only_resized_triangle(j_component)%line(i_component)) then
                    allocate(cells_target%neighbour_cells(j_component)%line(i_component)%cells, &
                        source=cells_source%neighbour_cells(j_component)%line(i_component)%cells)
                end if
            end do
        end do

        call visitable_cells_memento%save(cells_target%visitable_cells, cells_source%&
            visitable_cells)
    end subroutine save

    subroutine restore(cells_target, only_resized_triangle, visitable_cells_memento, cells_source)
        type(Cells_Wrapper), intent(inout) :: cells_target
        type(Logical_Line), intent(in) ::only_resized_triangle(:)
        class(Abstract_Visitable_Cells_Memento), intent(in) :: visitable_cells_memento
        type(Cells_Wrapper), intent(inout) :: cells_source

        integer :: i_component, j_component

        do j_component = 1, size(cells_target%neighbour_cells)
            do i_component = 1, size(cells_target%neighbour_cells(j_component)%line)
                if (only_resized_triangle(j_component)%line(i_component)) then
                    call cells_target%neighbour_cells(j_component)%line(i_component)%cells%reset()
                else
                    call cells_destroy(cells_target%neighbour_cells(j_component)%line(i_component)%&
                        cells)
                    allocate(cells_target%neighbour_cells(j_component)%line(i_component)%cells, &
                        source=cells_source%neighbour_cells(j_component)%line(i_component)%cells)
                end if
            end do
        end do
        call cells_destroy(cells_source%neighbour_cells)

        call visitable_cells_memento%restore(cells_target%visitable_cells, cells_target%&
            neighbour_cells, only_resized_triangle, cells_source%visitable_cells)
        call cells_destroy(cells_source%visitable_cells)
    end subroutine restore

end module procedures_cells_memento
