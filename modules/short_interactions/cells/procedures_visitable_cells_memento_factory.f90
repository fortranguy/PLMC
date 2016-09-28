module procedures_visitable_cells_memento_factory

use procedures_errors, only: error_exit
use classes_visitable_list, only: Abstract_Visitable_List
use classes_visitable_cells_memento, only: Abstract_Visitable_Cells_Memento, &
    Concrete_Visitable_Lists_Memento, Concrete_Visitable_Arrays_Memento, &
    Null_Visitable_Cells_Memento
use procedures_short_interactions_inquirers, only: list_is_linked_list, list_is_array

implicit none

private
public :: create, destroy

contains

    subroutine create(cells_memento, list_mold, needed)
        class(Abstract_Visitable_Cells_Memento), allocatable, intent(out) :: cells_memento
        class(Abstract_Visitable_List), intent(in) :: list_mold
        logical, intent(in) :: needed

        if (needed) then
            if (list_is_linked_list(list_mold)) then
                allocate(Concrete_Visitable_Lists_Memento :: cells_memento)
            else if (list_is_array(list_mold)) then
                allocate(Concrete_Visitable_Arrays_Memento :: cells_memento)
            else
                call error_exit("procedures_visitable_cells_memento_factory: create: "//&
                    "list_mold type is unknown")
            end if
        else
            allocate(Null_Visitable_Cells_Memento :: cells_memento)
        end if
    end subroutine create

    subroutine destroy(cells_memento)
        class(Abstract_Visitable_Cells_Memento), allocatable, intent(inout) :: cells_memento

        if (allocated(cells_memento)) deallocate(cells_memento)
    end subroutine destroy

end module procedures_visitable_cells_memento_factory
