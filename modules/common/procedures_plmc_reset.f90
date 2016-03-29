module procedures_plmc_reset

use class_mixture_total_moment, only: Abstract_Mixture_Total_Moment
use types_dipolar_interactions_wrapper, only: DES_Real_Pairs_Wrapper, Dipolar_Interactions_Wrapper

implicit none

private
public :: plmc_reset

interface plmc_reset
    module procedure :: reset_total_moments
    module procedure :: reset_dipolar
    module procedure :: reset_dipolar_real
end interface plmc_reset

contains

    subroutine reset_total_moments(total_moment)
        class(Abstract_Mixture_Total_Moment), intent(inout) :: total_moment

        call total_moment%reset()
    end subroutine reset_total_moments

    !> Some quantities may need to be reset to reflect the current configuration.
    subroutine reset_dipolar(dipolar_interactions)
        type(Dipolar_Interactions_Wrapper), intent(inout) :: dipolar_interactions

        call plmc_reset(dipolar_interactions%real_pairs)
        call dipolar_interactions%reci_weight%reset()
        call dipolar_interactions%reci_structure%reset()
        call dipolar_interactions%dlc_weight%reset()
        call dipolar_interactions%dlc_structures%reset()
    end subroutine reset_dipolar

    subroutine reset_dipolar_real(real_pairs)
        type(DES_Real_Pairs_Wrapper), intent(inout) :: real_pairs(:)

        integer :: i_component, j_component

        do j_component = 1, size(real_pairs)
            do i_component = 1, size(real_pairs(j_component)%line)
                call real_pairs(j_component)%line(i_component)%potential%reset()
            end do
        end do
    end subroutine reset_dipolar_real

end module procedures_plmc_reset
