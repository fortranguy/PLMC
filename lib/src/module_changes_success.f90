module module_changes_success

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_components

implicit none

private
public :: Concrete_Mixture_Changes_Counters_reset, Concrete_Mixture_Changes_Success_set

    type, public :: Concrete_Change_Counters
        integer :: num_hits
        integer :: num_success
    end type Concrete_Change_Counters

    type, public :: Concrete_Mixture_Changes_Counters
        type(Concrete_Change_Counters) :: moves(num_components)
        type(Concrete_Change_Counters) :: rotations(num_components)
        type(Concrete_Change_Counters) :: exchanges(num_components)
    end type Concrete_Mixture_Changes_Counters

    type, public :: Concrete_Changes_Success_Ratio
        real(DP) :: move
        real(DP) :: rotation
        real(DP) :: exchange
    end type Concrete_Changes_Success_Ratio

    type, public :: Concrete_Mixture_Changes_Success
        type(Concrete_Changes_Success_Ratio) :: ratios(num_components)
    end type Concrete_Mixture_Changes_Success

contains

    subroutine Concrete_Mixture_Changes_Counters_reset(changes_counters)
        type(Concrete_Mixture_Changes_Counters), intent(out) :: changes_counters

        integer :: i_component
        do i_component = 1, num_components
            call Concrete_Change_Counters_reset(changes_counters%moves(i_component))
            call Concrete_Change_Counters_reset(changes_counters%rotations(i_component))
            call Concrete_Change_Counters_reset(changes_counters%exchanges(i_component))
        end do
    end subroutine Concrete_Mixture_Changes_Counters_reset

    subroutine Concrete_Change_Counters_reset(change_counters)
        type(Concrete_Change_Counters), intent(out) :: change_counters

        change_counters%num_hits = 0
        change_counters%num_success = 0
    end subroutine Concrete_Change_Counters_reset

    subroutine Concrete_Mixture_Changes_Success_set(changes_success, changes_counters)
        type(Concrete_Mixture_Changes_Success), intent(out) :: changes_success
        type(Concrete_Mixture_Changes_Counters), intent(in) :: changes_counters

        integer :: i_component
        do i_component = 1, num_components
            changes_success%ratios(i_component)%move = &
                real(changes_counters%moves(i_component)%num_success, DP) / &
                real(changes_counters%moves(i_component)%num_hits, DP)
            changes_success%ratios(i_component)%rotation = &
                real(changes_counters%rotations(i_component)%num_success, DP) / &
                real(changes_counters%rotations(i_component)%num_hits, DP)
            changes_success%ratios(i_component)%exchange = &
                real(changes_counters%exchanges(i_component)%num_success, DP) / &
                real(changes_counters%exchanges(i_component)%num_hits, DP)
        end do
    end subroutine Concrete_Mixture_Changes_Success_set

end module module_changes_success
