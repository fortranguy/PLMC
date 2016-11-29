module procedures_dipolar_interactions_resetter

use classes_box_size_memento, only: Abstract_Box_Size_Memento
use classes_des_real_pair, only: Abstract_DES_Real_Pair
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper

implicit none

private
public :: reset, reset_real

contains

    !> Some dipolar accumulators may need to be reset to reflect the current configuration.
    subroutine reset(dipolar_interactions_static, reset_real_pair)
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: dipolar_interactions_static
        logical, intent(in) :: reset_real_pair

        call reset_real(dipolar_interactions_static%box_size_memento_real, &
            dipolar_interactions_static%real_pair, reset_real_pair)
        call dipolar_interactions_static%box_size_memento_reci%save()
        call dipolar_interactions_static%reci_weight%reset()
        call dipolar_interactions_static%reci_structure%reset()
        call dipolar_interactions_static%dlc_weight%reset()
        call dipolar_interactions_static%dlc_structures%reset()
    end subroutine reset

    subroutine reset_real(box_size_memento, real_pair, reset_real_pair)
        class(Abstract_Box_Size_Memento), intent(inout) :: box_size_memento
        class(Abstract_DES_Real_Pair), intent(inout) :: real_pair
        logical, intent(in) :: reset_real_pair

        if (reset_real_pair) then
            call box_size_memento%save()
            call real_pair%reset()
        end if
    end subroutine reset_real

end module procedures_dipolar_interactions_resetter
