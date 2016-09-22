module procedures_dipolar_interactions_resetter

use classes_box_volume_memento, only: Abstract_Box_Volume_Memento
use classes_des_real_pair, only: Abstract_DES_Real_Pair
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper

implicit none

private
public :: reset, reset_real

contains

    !> Some dipolar accumulators may need to be reset to reflect the current configuration.
    subroutine reset(dipolar_interactions_static, real_pair_must_be_reset)
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: dipolar_interactions_static
        logical, intent(in) :: real_pair_must_be_reset

        call reset_real(dipolar_interactions_static%box_volume_memento_real, &
            dipolar_interactions_static%real_pair, real_pair_must_be_reset)
        call dipolar_interactions_static%box_volume_memento_reci%save()
        call dipolar_interactions_static%reci_weight%reset()
        call dipolar_interactions_static%reci_structure%reset()
        call dipolar_interactions_static%dlc_weight%reset()
        call dipolar_interactions_static%dlc_structures%reset()
    end subroutine reset

    subroutine reset_real(box_volume_memento, real_pair, real_pair_must_be_reset)
        class(Abstract_Box_Volume_Memento), intent(inout) :: box_volume_memento
        class(Abstract_DES_Real_Pair), intent(inout) :: real_pair
        logical, intent(in) :: real_pair_must_be_reset

        if (real_pair_must_be_reset) then
            call box_volume_memento%save()
            call real_pair%reset()
        else
            !call real_pair%set_domain_max()
        end if
    end subroutine reset_real

end module procedures_dipolar_interactions_resetter
