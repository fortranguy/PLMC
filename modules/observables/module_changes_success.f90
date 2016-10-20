module module_changes_success

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_real_line, only: Real_Line

implicit none

private
public :: reset_counters, set_successes

interface reset_counters
    module procedure :: changes_counter_reset
    module procedure :: change_counters_reset
    module procedure :: change_counter_reset
end interface reset_counters

interface set_successes
    module procedure :: changes_success_set
    module procedure :: change_successes_set
    module procedure :: change_success_set
end interface set_successes

    type, public :: Concrete_Change_Counter
        integer :: num_hits = 0
        integer :: num_successes = 0
    end type Concrete_Change_Counter

    type, public :: Concrete_Changes_Counter
        type(Concrete_Change_Counter) :: translation
        type(Concrete_Change_Counter) :: rotation
        type(Concrete_Change_Counter) :: add
        type(Concrete_Change_Counter) :: remove
    end type Concrete_Changes_Counter

    type, public :: Concrete_Changes_Success
        real(DP) :: translation = 0._DP
        real(DP) :: rotation = 0._DP
        real(DP) :: add = 0._DP
        real(DP) :: remove = 0._DP
    end type Concrete_Changes_Success

    type, public :: Concrete_Change_Counter_Line
        type(Concrete_Change_Counter), allocatable :: line(:)
    end type Concrete_Change_Counter_Line

contains

    elemental subroutine changes_counter_reset(counter)
        type(Concrete_Changes_Counter), intent(inout) :: counter

        call reset_counters(counter%translation)
        call reset_counters(counter%rotation)
        call reset_counters(counter%add)
        call reset_counters(counter%remove)
    end subroutine changes_counter_reset

    elemental subroutine change_counters_reset(counters)
        type(Concrete_Change_Counter_Line), intent(inout) :: counters

        call reset_counters(counters%line)
    end subroutine change_counters_reset

    elemental subroutine change_counter_reset(counter)
        type(Concrete_Change_Counter), intent(inout) :: counter

        counter%num_hits = 0
        counter%num_successes = 0
    end subroutine change_counter_reset

    elemental subroutine changes_success_set(success, counter)
        type(Concrete_Changes_Success), intent(inout) :: success
        type(Concrete_Changes_Counter), intent(in) :: counter

        success%translation = calculate_ratio(counter%translation)
        success%rotation = calculate_ratio(counter%rotation)
        success%add = calculate_ratio(counter%add)
        success%remove = calculate_ratio(counter%remove)
    end subroutine changes_success_set

    elemental subroutine change_successes_set(successes, counters)
        type(Real_Line), intent(inout) :: successes
        type(Concrete_Change_Counter_Line), intent(in) :: counters

        successes%line = calculate_ratio(counters%line)
    end subroutine change_successes_set

    elemental subroutine change_success_set(success, counter)
        real(DP), intent(inout) :: success
        type(Concrete_Change_Counter), intent(in) :: counter

        success = calculate_ratio(counter)
    end subroutine change_success_set

    elemental real(DP) function calculate_ratio(counter) result(ratio)
        type(Concrete_Change_Counter), intent(in) :: counter

        if (counter%num_hits > 0) then
            ratio = real(counter%num_successes, DP) / real(counter%num_hits, DP)
        else
            ratio = 0._DP
        end if
    end function calculate_ratio

end module module_changes_success
