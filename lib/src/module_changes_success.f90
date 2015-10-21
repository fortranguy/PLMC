module module_changes_success

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public :: Concrete_Changes_Counter_reset, Concrete_Changes_Counter_set

    type, public :: Concrete_Change_Counter
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

contains

    pure subroutine Concrete_Changes_Counter_reset(changes_counter)
        type(Concrete_Changes_Counter), intent(out) :: changes_counter

        call reset(changes_counter%move)
        call reset(changes_counter%rotation)
        call reset(changes_counter%exchange)
    end subroutine Concrete_Changes_Counter_reset

    pure subroutine reset(changes_counter)
        type(Concrete_Change_Counter), intent(out) :: changes_counter

        changes_counter%num_hits = 0
        changes_counter%num_success = 0
    end subroutine reset

    pure subroutine Concrete_Changes_Counter_set(changes_success, changes_counter)
        type(Concrete_Changes_Success), intent(out) :: changes_success
        type(Concrete_Changes_Counter), intent(in) :: changes_counter

        changes_success%move = get_ratio(changes_counter%move)
        changes_success%rotation = get_ratio(changes_counter%rotation)
        changes_success%exchange = get_ratio(changes_counter%exchange)
    end subroutine Concrete_Changes_Counter_set

    pure real(DP) function get_ratio(change_counter)
        type(Concrete_Change_Counter), intent(in) :: change_counter

        if (change_counter%num_hits > 0) then
            get_ratio = real(change_counter%num_success, DP) / &
                real(change_counter%num_hits, DP)
        else
            get_ratio = 0._DP
        end if
    end function get_ratio

end module module_changes_success
