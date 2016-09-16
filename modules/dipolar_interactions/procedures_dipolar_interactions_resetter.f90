module procedures_dipolar_interactions_resetter

use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper

implicit none

private
public :: reset

contains

    !> Some dipolar accumulators may need to be reset to reflect the current configuration.
    subroutine reset(dipolar_interactions_static, real_pair_must_be_reset)
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: dipolar_interactions_static
        logical, intent(in) :: real_pair_must_be_reset

        if (real_pair_must_be_reset) then
            call dipolar_interactions_static%box_volume_memento_real%save()
            call dipolar_interactions_static%real_pair%reset()
        end if
        call dipolar_interactions_static%box_volume_memento_reci%save()
        call dipolar_interactions_static%reci_weight%reset()
        call dipolar_interactions_static%reci_structure%reset()
        call dipolar_interactions_static%dlc_weight%reset()
        call dipolar_interactions_static%dlc_structures%reset()
    end subroutine reset

end module procedures_dipolar_interactions_resetter
