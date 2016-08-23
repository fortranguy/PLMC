module classes_visitable_cells_memento

use types_neighbour_cells_wrapper, only: Neighbour_Cells_Line
use classes_visitable_cells, only: Abstract_Visitable_Cells, Null_Visitable_Cells
use procedures_visitable_cells_factory, only: visitable_cells_destroy => destroy

implicit none

private

    type, abstract, public :: Abstract_Visitable_Cells_Memento
    contains
        procedure(Abstract_save), deferred, nopass :: save
        procedure(Abstract_restore), deferred, nopass :: restore
    end type Abstract_Visitable_Cells_Memento

    abstract interface

        subroutine Abstract_save(visitable_cells_target, visitable_cells_source)
        import :: Abstract_Visitable_Cells
            class(Abstract_Visitable_Cells), allocatable, intent(out) :: &
                visitable_cells_target(:, :)
            class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)
        end subroutine Abstract_save

        subroutine Abstract_restore(visitable_cells_target, visitable_cells_source, neighbour_cells)
        import :: Neighbour_Cells_Line, Abstract_Visitable_Cells
            class(Abstract_Visitable_Cells), allocatable, intent(inout) :: &
                visitable_cells_target(:, :)
            class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)
            type(Neighbour_Cells_Line), intent(in) :: neighbour_cells(:)
        end subroutine Abstract_restore

    end interface

    type, extends(Abstract_Visitable_Cells_Memento), public :: Concrete_Visitable_Lists_Memento
    contains
        procedure, nopass :: save => Lists_save
        procedure, nopass :: restore => Lists_restore
    end type Concrete_Visitable_Lists_Memento

    type, extends(Abstract_Visitable_Cells_Memento), public :: Concrete_Visitable_Arrays_Memento
    contains
        procedure, nopass :: save => Arrays_save
        procedure, nopass :: restore => Arrays_restore
    end type Concrete_Visitable_Arrays_Memento

    type, extends(Abstract_Visitable_Cells_Memento), public :: Null_Visitable_Cells_Memento
    contains
        procedure, nopass :: save => Null_save
        procedure, nopass :: restore => Null_restore
    end type Null_Visitable_Cells_Memento

contains

!implementation Concrete_Visitable_Lists_Memento

    subroutine Lists_save(visitable_cells_target, visitable_cells_source)
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: visitable_cells_target(:, :)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)

        allocate(Null_Visitable_Cells :: visitable_cells_target(size(visitable_cells_source, 1), &
            size(visitable_cells_source, 2)))
    end subroutine Lists_save

    !> @note Instead of copying linked-lists (i.e. array of [[Abstract_Visitable_List]]),
    !> [[Lists_restore]] rebuilds them.
    subroutine Lists_restore(visitable_cells_target, visitable_cells_source, neighbour_cells)
        class(Abstract_Visitable_Cells), allocatable, intent(inout) :: visitable_cells_target(:, :)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)
        type(Neighbour_Cells_Line), intent(in) :: neighbour_cells(:)

        integer :: i_component, j_component
        integer :: i_pair, j_pair

        do j_component = 1, size(visitable_cells_source, 2)
            do i_component = 1, size(visitable_cells_source, 1)
                j_pair = maxval([i_component, j_component])
                i_pair = minval([i_component, j_component])
                call visitable_cells_target(i_component, j_component)%set(neighbour_cells(j_pair)%&
                    line(i_pair)%cells)
                call visitable_cells_target(i_component, j_component)%reset()
            end do
        end do
    end subroutine Lists_restore

!end implementation Concrete_Visitable_Lists_Memento

!implementation Concrete_Visitable_Arrays_Memento

    subroutine Arrays_save(visitable_cells_target, visitable_cells_source)
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: visitable_cells_target(:, :)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)

        allocate(visitable_cells_target, source=visitable_cells_source)
    end subroutine Arrays_save

    subroutine Arrays_restore(visitable_cells_target, visitable_cells_source, neighbour_cells)
        class(Abstract_Visitable_Cells), allocatable, intent(inout) :: visitable_cells_target(:, :)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)
        type(Neighbour_Cells_Line), intent(in) :: neighbour_cells(:)

        integer :: i_component, j_component
        integer :: i_pair, j_pair

        call visitable_cells_destroy(visitable_cells_target)
        allocate(visitable_cells_target, source=visitable_cells_source)
        do j_component = 1, size(visitable_cells_target, 2)
            do i_component = 1, size(visitable_cells_target, 1)
                j_pair = maxval([i_component, j_component])
                i_pair = minval([i_component, j_component])
                call visitable_cells_target(i_component, j_component)%set(neighbour_cells(j_pair)%&
                    line(i_pair)%cells)
            end do
        end do
    end subroutine Arrays_restore

!end implementation Concrete_Visitable_Arrays_Memento

!implementation Null_Visitable_Cells_Memento

    subroutine Null_save(visitable_cells_target, visitable_cells_source)
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: visitable_cells_target(:, :)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)
    end subroutine Null_save

    subroutine Null_restore(visitable_cells_target, visitable_cells_source, neighbour_cells)
        class(Abstract_Visitable_Cells), allocatable, intent(inout) :: visitable_cells_target(:, :)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells_source(:, :)
        type(Neighbour_Cells_Line), intent(in) :: neighbour_cells(:)
    end subroutine Null_restore

!end implementation Null_Visitable_Cells_Memento

end module classes_visitable_cells_memento
