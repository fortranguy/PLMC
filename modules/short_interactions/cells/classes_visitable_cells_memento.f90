module classes_visitable_cells_memento

use types_logical_line, only: Logical_Line
use classes_neighbour_cells, only: Neighbour_Cells_Line
use classes_visitable_cells, only: Abstract_Visitable_Cells, Null_Visitable_Cells
use procedures_visitable_cells_factory, only: visitable_cells_destroy => destroy

implicit none

private

    !> @bug [[classes_visitable_cells_memento:Abstract_save]] must be nopass: gfortran bug?
    type, abstract, public :: Abstract_Visitable_Cells_Memento
    contains
        procedure(Abstract_save), deferred :: save
        procedure(Abstract_restore), deferred, nopass :: restore
    end type Abstract_Visitable_Cells_Memento

    abstract interface

        subroutine Abstract_save(this, visitable_cells_target, visitable_cells_source)
        import :: Abstract_Visitable_Cells, Abstract_Visitable_Cells_Memento
            class(Abstract_Visitable_Cells_Memento), intent(in) :: this
            class(Abstract_Visitable_Cells), allocatable, intent(out) :: &
                visitable_cells_target(:, :)
            class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)
        end subroutine Abstract_save

        subroutine Abstract_restore(visitable_cells_target, neighbour_cells, only_resized_triangle,&
            visitable_cells_source)
        import :: Logical_Line, Neighbour_Cells_Line, Abstract_Visitable_Cells
            class(Abstract_Visitable_Cells), allocatable, intent(inout) :: &
                visitable_cells_target(:, :)
            type(Neighbour_Cells_Line), intent(in) :: neighbour_cells(:)
            type(Logical_Line), intent(in) ::only_resized_triangle(:)
            class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)
        end subroutine Abstract_restore

    end interface

    type, extends(Abstract_Visitable_Cells_Memento), public :: Concrete_Visitable_Lists_Memento
    contains
        procedure :: save => Lists_save
        procedure, nopass :: restore => Lists_restore
    end type Concrete_Visitable_Lists_Memento

    type, extends(Abstract_Visitable_Cells_Memento), public :: Concrete_Visitable_Arrays_Memento
    contains
        procedure :: save => Arrays_save
        procedure, nopass :: restore => Arrays_restore
    end type Concrete_Visitable_Arrays_Memento

    type, extends(Abstract_Visitable_Cells_Memento), public :: Null_Visitable_Cells_Memento
    contains
        procedure :: save => Null_save
        procedure, nopass :: restore => Null_restore
    end type Null_Visitable_Cells_Memento

contains

!implementation Concrete_Visitable_Lists_Memento

    subroutine Lists_save(this, visitable_cells_target, visitable_cells_source)
        class(Concrete_Visitable_Lists_Memento), intent(in) :: this
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: visitable_cells_target(:, :)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)

        allocate(Null_Visitable_Cells :: visitable_cells_target(size(visitable_cells_source, 1), &
            size(visitable_cells_source, 2)))
    end subroutine Lists_save

    !> @note Instead of copying linked-lists (i.e. array of [[Abstract_Visitable_List]]),
    !> [[classes_visitable_cells_memento:Lists_restore]] rebuilds them.
    subroutine Lists_restore(visitable_cells_target, neighbour_cells, only_resized_triangle, &
        visitable_cells_source)
        class(Abstract_Visitable_Cells), allocatable, intent(inout) :: visitable_cells_target(:, :)
        type(Neighbour_Cells_Line), intent(in) :: neighbour_cells(:)
        type(Logical_Line), intent(in) ::only_resized_triangle(:)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)

        integer :: i_component, j_component
        integer :: i_pair, j_pair

        do j_component = 1, size(visitable_cells_source, 2)
            do i_component = 1, size(visitable_cells_source, 1)
                j_pair = maxval([i_component, j_component])
                i_pair = minval([i_component, j_component])
                call visitable_cells_target(i_component, j_component)%&
                    target(neighbour_cells(j_pair)%line(i_pair)%cells)
                if (only_resized_triangle(j_pair)%line(i_pair)) cycle
                call visitable_cells_target(i_component, j_component)%reset()
            end do
        end do
    end subroutine Lists_restore

!end implementation Concrete_Visitable_Lists_Memento

!implementation Concrete_Visitable_Arrays_Memento

    subroutine Arrays_save(this, visitable_cells_target, visitable_cells_source)
        class(Concrete_Visitable_Arrays_Memento), intent(in) :: this
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: visitable_cells_target(:, :)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)

        allocate(visitable_cells_target, source=visitable_cells_source)
    end subroutine Arrays_save

    subroutine Arrays_restore(visitable_cells_target, neighbour_cells, only_resized_triangle, &
        visitable_cells_source)
        class(Abstract_Visitable_Cells), allocatable, intent(inout) :: visitable_cells_target(:, :)
        type(Neighbour_Cells_Line), intent(in) :: neighbour_cells(:)
        type(Logical_Line), intent(in) ::only_resized_triangle(:) !useful?
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)

        integer :: i_component, j_component
        integer :: i_pair, j_pair
        logical :: only_resized

        only_resized = .true.
        do j_component = 1, size(only_resized_triangle)
            only_resized = only_resized .and. all(only_resized_triangle(j_component)%line)
            if (.not.only_resized) exit
        end do
        if (.not.only_resized) then
            call visitable_cells_destroy(visitable_cells_target)
            allocate(visitable_cells_target, source=visitable_cells_source)
        end if
        do j_component = 1, size(visitable_cells_target, 2)
            do i_component = 1, size(visitable_cells_target, 1)
                j_pair = maxval([i_component, j_component])
                i_pair = minval([i_component, j_component])
                call visitable_cells_target(i_component, j_component)%&
                    target(neighbour_cells(j_pair)%line(i_pair)%cells)
            end do
        end do
    end subroutine Arrays_restore

!end implementation Concrete_Visitable_Arrays_Memento

!implementation Null_Visitable_Cells_Memento

    subroutine Null_save(this, visitable_cells_target, visitable_cells_source)
        class(Null_Visitable_Cells_Memento), intent(in) :: this
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: visitable_cells_target(:, :)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)
    end subroutine Null_save

    subroutine Null_restore(visitable_cells_target, neighbour_cells, only_resized_triangle, &
        visitable_cells_source)
        class(Abstract_Visitable_Cells), allocatable, intent(inout) :: visitable_cells_target(:, :)
        type(Neighbour_Cells_Line), intent(in) :: neighbour_cells(:)
        type(Logical_Line), intent(in) ::only_resized_triangle(:)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)
    end subroutine Null_restore

!end implementation Null_Visitable_Cells_Memento

end module classes_visitable_cells_memento
