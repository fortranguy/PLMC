module procedures_dipolar_neighbourhoods_visitors_factory

use types_physical_model_wrapper, only: Physical_Model_Wrapper
use classes_dipolar_neighbourhoods_visitor, only: Abstract_Dipolar_Neighbourhoods_Visitor, &
    Concrete_Dipolar_Neighbourhoods_Visitor, Null_Dipolar_Neighbourhoods_Visitor

implicit none

private
public :: create, destroy

contains

    subroutine create(dipolar_neighbourhoods_visitors, physical_model, make_dipoles_graph)
        class(Abstract_Dipolar_Neighbourhoods_Visitor), allocatable, intent(out) :: &
            dipolar_neighbourhoods_visitors(:)
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        logical, intent(in) :: make_dipoles_graph

        integer :: i_box

        if (make_dipoles_graph) then
            allocate(Concrete_Dipolar_Neighbourhoods_Visitor :: &
                dipolar_neighbourhoods_visitors(size(physical_model%environment%periodic_boxes)))
        else
            allocate(Null_Dipolar_Neighbourhoods_Visitor :: &
                dipolar_neighbourhoods_visitors(size(physical_model%environment%periodic_boxes)))
        end if

        do i_box = 1, size(dipolar_neighbourhoods_visitors)
            call dipolar_neighbourhoods_visitors(i_box)%&
                construct(physical_model%mixture%components(:, i_box), physical_model%&
                short_interactions%cells(i_box)%visitable_cells)
        end do
    end subroutine create

    subroutine destroy(dipolar_neighbourhoods_visitors)
        class(Abstract_Dipolar_Neighbourhoods_Visitor), allocatable, intent(inout) :: &
            dipolar_neighbourhoods_visitors(:)

        integer :: i_box

        if (allocated(dipolar_neighbourhoods_visitors)) then
            do i_box = size(dipolar_neighbourhoods_visitors), 1, -1
                call dipolar_neighbourhoods_visitors(i_box)%destroy()
            end do
            deallocate(dipolar_neighbourhoods_visitors)
        end if
    end subroutine destroy

end module procedures_dipolar_neighbourhoods_visitors_factory
