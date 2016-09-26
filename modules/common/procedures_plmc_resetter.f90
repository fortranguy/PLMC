module procedures_plmc_resetter

use procedures_short_interactions_resetter, only: short_interactions_reset => reset
use procedures_dipolar_interactions_resetter, only: dipolar_interactions_reset => reset
use types_physical_model_wrapper, only: Physical_Model_Wrapper

implicit none

private
public :: plmc_reset

contains

    subroutine plmc_reset(physical_model)
        type(Physical_Model_Wrapper), intent(inout) :: physical_model

        logical :: reset_real_pair

        call physical_model%mixture%total_moment%reset()
        call short_interactions_reset(physical_model%short_interactions%neighbour_cells)
        call short_interactions_reset(physical_model%short_interactions%visitable_cells)
        reset_real_pair = .true.
        call dipolar_interactions_reset(physical_model%dipolar_interactions_static, reset_real_pair)
    end subroutine plmc_reset

end module procedures_plmc_resetter
