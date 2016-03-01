module module_changes_success

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_line_observables, only: Concrete_Line_Observables

implicit none

private
public :: changes_counter_reset, changes_counter_set, switches_counters_reset, switches_counters_set

    type :: Concrete_Change_Counter
        integer :: num_hits
        integer :: num_success
    end type Concrete_Change_Counter

    type, public :: Concrete_Changes_Counter
        type(Concrete_Change_Counter) :: move
        type(Concrete_Change_Counter) :: rotation
        type(Concrete_Change_Counter) :: exchange
    end type Concrete_Changes_Counter

    type, public :: Concrete_Changes_Success
        real(DP) :: move
        real(DP) :: rotation
        real(DP) :: exchange
    end type Concrete_Changes_Success

    type, public :: Concrete_Switch_Counters
        type(Concrete_Change_Counter), allocatable :: with_components(:)
    end type Concrete_Switch_Counters

contains

    pure elemental subroutine changes_counter_reset(changes_counter)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counter

        call reset(changes_counter%move)
        call reset(changes_counter%rotation)
        call reset(changes_counter%exchange)
    end subroutine changes_counter_reset

    pure elemental subroutine switches_counters_reset(switches_counters)
        type(Concrete_Switch_Counters), intent(inout) :: switches_counters

        call reset(switches_counters%with_components)
    end subroutine switches_counters_reset

    pure elemental subroutine reset(changes_counter)
        type(Concrete_Change_Counter), intent(inout) :: changes_counter

        changes_counter%num_hits = 0
        changes_counter%num_success = 0
    end subroutine reset

    pure elemental subroutine changes_counter_set(changes_success, changes_counter)
        type(Concrete_Changes_Success), intent(inout) :: changes_success
        type(Concrete_Changes_Counter), intent(in) :: changes_counter

        changes_success%move = get_ratio(changes_counter%move)
        changes_success%rotation = get_ratio(changes_counter%rotation)
        changes_success%exchange = get_ratio(changes_counter%exchange)
    end subroutine changes_counter_set

    pure elemental subroutine switches_counters_set(switches_successes, switches_counters)
        type(Concrete_Line_Observables), intent(inout) :: switches_successes
        type(Concrete_Switch_Counters), intent(in) :: switches_counters

        switches_successes%with_components = get_ratio(switches_counters%with_components)
    end subroutine switches_counters_set

    pure elemental real(DP) function get_ratio(change_counter)
        type(Concrete_Change_Counter), intent(in) :: change_counter

        if (change_counter%num_hits > 0) then
            get_ratio = real(change_counter%num_success, DP) / real(change_counter%num_hits, DP)
        else
            get_ratio = 0._DP
        end if
    end function get_ratio

end module module_changes_success
