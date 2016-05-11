module classes_move_tuner

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: warning_continue
use procedures_checks, only: check_positive, check_ratio
use classes_moved_coordinates, only: Abstract_Moved_Coordinates
use types_move_tuner_parameters, only: Concrete_Move_Tuner_Parameters

implicit none

private

    type, abstract, public :: Abstract_Move_Tuner
    private
        class(Abstract_Moved_Coordinates), pointer :: moved_coordinates => null()
        integer :: accumulation_period = 0
        real(DP) :: accumulated_success_ratio = 0._DP
        real(DP) :: success_ratio_min = 0._DP, success_ratio_max = 0._DP
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: tune => Abstract_tune
    end type Abstract_Move_Tuner

    type, extends(Abstract_Move_Tuner), public :: Concrete_Move_Tuner

    end type Concrete_Move_Tuner

    type, extends(Abstract_Move_Tuner), public :: Null_Move_Tuner
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: tune => Null_tune
    end type Null_Move_Tuner

contains

!implementation Abstract_Move_Tuner

    subroutine Abstract_construct(this, moved_coordinates, parameters, num_tuning_steps)
        class(Abstract_Move_Tuner), intent(out) :: this
        class(Abstract_Moved_Coordinates), target, intent(in) :: moved_coordinates
        type(Concrete_Move_Tuner_Parameters), intent(in) :: parameters
        integer, intent(in) :: num_tuning_steps

        this%moved_coordinates => moved_coordinates
        call check_positive("Abstract_Move_Tuner: construct", "parameters%accumulation_period", &
            parameters%accumulation_period)
        if (num_tuning_steps < parameters%accumulation_period) then
            call warning_continue("Abstract_Move_Tuner: construct: "//&
                "num_tuning_steps < accumulation_period")
        end if
        this%accumulation_period = parameters%accumulation_period
        this%accumulated_success_ratio = 0._DP
        call check_ratio("Abstract_Move_Tuner: construct", "parameters%wanted_success_ratio", &
            parameters%wanted_success_ratio)
        call check_positive("Abstract_Move_Tuner: construct", "parameters%tolerance", &
            parameters%tolerance)
        this%success_ratio_min = parameters%wanted_success_ratio - parameters%tolerance
        this%success_ratio_max = parameters%wanted_success_ratio + parameters%tolerance
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Move_Tuner), intent(inout) :: this

        this%moved_coordinates => null()
    end subroutine Abstract_destroy

    subroutine Abstract_tune(this, tuned, i_step, success_ratio)
        class(Abstract_Move_Tuner), intent(inout) :: this
        logical, intent(out) :: tuned
        integer, intent(in) :: i_step
        real(DP), intent(in) :: success_ratio

        real(DP) :: averaged_success_ratio

        tuned = .false.
        this%accumulated_success_ratio = this%accumulated_success_ratio + success_ratio
        if (mod(i_step, this%accumulation_period) == 0) then
            averaged_success_ratio = this%accumulated_success_ratio / &
                real(this%accumulation_period, DP) ! assuming i_step++;
            if (averaged_success_ratio < this%success_ratio_min) then
                call this%moved_coordinates%decrease_delta()
            else if (this%success_ratio_max < averaged_success_ratio) then
                call this%moved_coordinates%increase_delta()
            else
                tuned = .true. ! sufficient condition?
                return
            end if
            this%accumulated_success_ratio = 0._DP
        end if
    end subroutine Abstract_tune

!end implementation Abstract_Move_Tuner

!implementation Null_Move_Tuner

    subroutine Null_construct(this, moved_coordinates, parameters, num_tuning_steps)
        class(Null_Move_Tuner), intent(out) :: this
        class(Abstract_Moved_Coordinates), target, intent(in) :: moved_coordinates
        type(Concrete_Move_Tuner_Parameters), intent(in) :: parameters
        integer, intent(in) :: num_tuning_steps
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Move_Tuner), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_tune(this, tuned, i_step, success_ratio)
        class(Null_Move_Tuner), intent(inout) :: this
        logical, intent(out) :: tuned
        integer, intent(in) :: i_step
        real(DP), intent(in) :: success_ratio
        tuned = .true.
    end subroutine Null_tune

!end implementation Null_Move_Tuner

end module classes_move_tuner
