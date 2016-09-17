module procedures_cells_memento

use types_logical_line, only: Concrete_Logical_Line
use types_neighbour_cells_wrapper, only: Neighbour_Cells_Line
use classes_visitable_cells, only: Abstract_Visitable_Cells
use procedures_cells_factory, only: cells_destroy => destroy
use classes_visitable_cells_memento, only: Abstract_Visitable_Cells_Memento

implicit none

private
public :: save, restore

contains

    subroutine save(visitable_cells_memento, neighbour_cells_target, visitable_cells_target, &
        neighbour_cells_source, visitable_cells_source)
        class(Abstract_Visitable_Cells_Memento), intent(in) :: visitable_cells_memento
        type(Neighbour_Cells_Line), allocatable, intent(out) :: neighbour_cells_target(:)
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: visitable_cells_target(:, :)
        type(Neighbour_Cells_Line), allocatable, intent(in) :: neighbour_cells_source(:)
        class(Abstract_Visitable_Cells), allocatable, intent(in) :: visitable_cells_source(:, :)

        integer :: i_component, j_component

        if (allocated(neighbour_cells_target)) deallocate(neighbour_cells_target)
        allocate(neighbour_cells_target(size(neighbour_cells_source)))
        do j_component = 1, size(neighbour_cells_target)
            allocate(neighbour_cells_target(j_component)%line(j_component))
            do i_component = 1, size(neighbour_cells_target(j_component)%line)
                allocate(neighbour_cells_target(j_component)%line(i_component)%cells, &
                    source=neighbour_cells_source(j_component)%line(i_component)%cells)
            end do
        end do
        if (allocated(visitable_cells_target)) deallocate(visitable_cells_target)
        call visitable_cells_memento%save(visitable_cells_target, visitable_cells_source)
    end subroutine save

    subroutine restore(visitable_cells_memento, neighbour_cells_target, visitable_cells_target, &
        neighbour_cells_source, only_resized_triangle, visitable_cells_source)
        class(Abstract_Visitable_Cells_Memento), intent(in) :: visitable_cells_memento
        type(Neighbour_Cells_Line), allocatable, intent(out) :: neighbour_cells_target(:)
        class(Abstract_Visitable_Cells), allocatable, intent(inout) :: visitable_cells_target(:, :)
        type(Neighbour_Cells_Line), intent(in) :: neighbour_cells_source(:)
        type(Concrete_Logical_Line), intent(in) ::only_resized_triangle(:)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)

        integer :: i_component, j_component

        call cells_destroy(neighbour_cells_target)
        allocate(neighbour_cells_target(size(neighbour_cells_source)))
        do j_component = 1, size(neighbour_cells_source)
            allocate(neighbour_cells_target(j_component)%line(j_component))
            do i_component = 1, size(neighbour_cells_source(j_component)%line)
                allocate(neighbour_cells_target(j_component)%line(i_component)%&
                    cells, source=neighbour_cells_source(j_component)%line(i_component)%cells)
            end do
        end do
        call visitable_cells_memento%restore(visitable_cells_target, neighbour_cells_target,&
            only_resized_triangle, visitable_cells_source)
    end subroutine restore

end module procedures_cells_memento
