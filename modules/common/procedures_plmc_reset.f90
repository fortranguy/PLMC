module procedures_plmc_reset

use types_long_interactions_wrapper, only: Ewald_Real_Pairs_Wrapper, Ewald_Reci_Component_Wrapper, &
    Long_Interactions_Wrapper

implicit none

private
public :: plmc_reset

interface plmc_reset
    module procedure :: reset_long
    module procedure :: reset_long_real
    module procedure :: reset_long_reci
end interface plmc_reset

contains

    !> When using Ewald, some quantities may need to be reset to reflect the actual configuration.
    subroutine reset_long(long_interactions)
        type(Long_Interactions_Wrapper), intent(inout) :: long_interactions

        call plmc_reset(long_interactions%real_pairs)
        call long_interactions%reci_weight%reset()
        call plmc_reset(long_interactions%reci_components)
    end subroutine reset_long

    !> To ensure volume changes were taken into account
    subroutine reset_long_real(real_pairs)
        type(Ewald_Real_Pairs_Wrapper), intent(inout) :: real_pairs(:)

        integer :: j_component, i_component

        do j_component = 1, size(real_pairs)
            do i_component = 1, size(real_pairs(j_component)%with_components)
                call real_pairs(j_component)%with_components(i_component)%real_pair%reset()
            end do
        end do
    end subroutine reset_long_real

    !> To ensure components changes were taken into account
    subroutine reset_long_reci(reci_components)
        type(Ewald_Reci_Component_Wrapper), intent(inout) :: reci_components(:)

        integer :: i_component

        do i_component = 1, size(reci_components)
            call reci_components(i_component)%reci_structure%reset()
        end do
    end subroutine reset_long_reci

end module procedures_plmc_reset
