module procedures_plmc_resetter

use procedures_short_interactions_resetter, only: short_interactions_reset => reset
use procedures_dipolar_interactions_resetter, only: dipolar_interactions_reset => reset
use types_physical_model_wrapper, only: Physical_Model_Wrapper

implicit none

private
public :: plmc_reset

contains

    !> @note average_num_particles%set() will be useless at the end of the run.
    subroutine plmc_reset(physical_model)
        type(Physical_Model_Wrapper), intent(inout) :: physical_model

        integer :: i_box, i_component
        logical :: reset_real_pair

        do i_box = 1, size(physical_model%mixture%gemc_components, 2)
            do i_component = 1, size(physical_model%mixture%gemc_components, 1)
                call physical_model%mixture%gemc_components(i_component, i_box)%average_num_particles%&
                    set()
            end do
        end do
        call physical_model%mixture%total_moment%reset()
        call short_interactions_reset(physical_model%short_interactions%neighbour_cells)
        call short_interactions_reset(physical_model%short_interactions%visitable_cells)
        reset_real_pair = .true.
        call dipolar_interactions_reset(physical_model%dipolar_interactions_static, reset_real_pair)
    end subroutine plmc_reset

end module procedures_plmc_resetter
