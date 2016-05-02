module procedures_plmc_reset

use types_neighbour_cells_wrapper, only: Neighbour_Cells_Line
use classes_visitable_cells, only: Abstract_Visitable_Cells
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use types_physical_model_wrapper, only: Physical_Model_Wrapper

implicit none

private
public :: plmc_reset

contains

    subroutine plmc_reset(physical_model)
        type(Physical_Model_Wrapper), intent(inout) :: physical_model

        call physical_model%mixture%total_moment%reset()
        call reset_neighbour_cells(physical_model%short_interactions%neighbour_cells)
        call reset_visitable_cells(physical_model%short_interactions%visitable_cells)
        call reset_dipolar(physical_model%dipolar_interactions)
    end subroutine plmc_reset

    subroutine reset_neighbour_cells(cells)
        type(Neighbour_Cells_Line), intent(inout) :: cells(:)

        integer :: i_component, j_component

        do j_component = 1, size(cells)
            do i_component = 1, size(cells(j_component)%line)
                call cells(j_component)%line(i_component)%cells%reset()
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

    !> Some dipolar accumulators may need to be reset to reflect the current configuration.
    subroutine reset_dipolar(dipolar_interactions)
        type(Dipolar_Interactions_Wrapper), intent(inout) :: dipolar_interactions

        call dipolar_interactions%real_pair%reset()
        call dipolar_interactions%reci_weight%reset()
        call dipolar_interactions%reci_structure%reset()
        call dipolar_interactions%dlc_weight%reset()
        call dipolar_interactions%dlc_structures%reset()
    end subroutine reset_dipolar

end module procedures_plmc_reset
