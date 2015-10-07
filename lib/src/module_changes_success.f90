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

        call reset(changes_counters%moves)
        call reset(changes_counters%rotations)
        call reset(changes_counters%exchanges)
    end subroutine Concrete_Mixture_Changes_Counters_reset

    elemental subroutine reset(change_counters)
        type(Concrete_Change_Counters), intent(out) :: change_counters

        change_counters%num_hits = 0
        change_counters%num_success = 0
    end subroutine reset

    subroutine Concrete_Mixture_Changes_Success_set(changes_success, changes_counters)
        type(Concrete_Mixture_Changes_Success), intent(out) :: changes_success
        type(Concrete_Mixture_Changes_Counters), intent(in) :: changes_counters

        changes_success%ratios(:)%move = get_ratio(changes_counters%moves)
        changes_success%ratios(:)%rotation = get_ratio(changes_counters%rotations)
        changes_success%ratios(:)%exchange = get_ratio(changes_counters%exchanges)
    end subroutine Concrete_Mixture_Changes_Success_set

    elemental real(DP) function get_ratio(change_counters)
        type(Concrete_Change_Counters), intent(in) :: change_counters

        if (change_counters%num_hits > 0) then
            get_ratio = real(change_counters%num_success, DP) / &
                real(change_counters%num_hits, DP)
        else
            get_ratio = 0._DP
        end if
    end function get_ratio

end module module_changes_success
