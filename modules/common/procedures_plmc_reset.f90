module procedures_plmc_reset

use types_long_interactions_wrapper, only: Ewald_Real_Pairs_Wrapper, Long_Interactions_Wrapper

implicit none

private
public :: plmc_reset

interface plmc_reset
    module procedure :: reset_long
    module procedure :: reset_long_real
end interface plmc_reset

contains

    !> When using Ewald, some quantities may need to be reset to reflect the actual configuration.
    subroutine reset_long(long_interactions)
        type(Long_Interactions_Wrapper), intent(inout) :: long_interactions

        call plmc_reset(long_interactions%real_pairs)
        call long_interactions%reci_weight%reset()
        call long_interactions%reci_structure%reset()
    end subroutine reset_long

    subroutine reset_long_real(real_pairs)
        type(Ewald_Real_Pairs_Wrapper), intent(inout) :: real_pairs(:)

        integer :: j_component, i_component

        do j_component = 1, size(real_pairs)
            do i_component = 1, size(real_pairs(j_component)%with_components)
                call real_pairs(j_component)%with_components(i_component)%real_pair%reset()
            end do
        end do
    end subroutine reset_long_real

end module procedures_plmc_reset
