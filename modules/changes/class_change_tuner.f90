module class_change_tuner

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: warning_continue
use procedures_checks, only: check_positive, check_ratio
use module_plmc_iterations, only: num_tuning_steps
use class_changed_coordinates, only: Abstract_Changed_Coordinates

implicit none

private

    type, public :: Concrete_Change_Tuner_Parameters
        integer :: accumulation_period
        real(DP) :: wanted_success_ratio
        real(DP) :: tolerance
    end type Concrete_Change_Tuner_Parameters

    type, abstract, public :: Abstract_Change_Tuner
    private
        class(Abstract_Changed_Coordinates), pointer :: changed_coordinates => null()
        integer :: accumulation_period
        real(DP) :: accumulated_success_ratio
        real(DP) :: success_ratio_min, success_ratio_max
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: tune => Abstract_tune
    end type Abstract_Change_Tuner

    type, extends(Abstract_Change_Tuner), public :: Concrete_Change_Tuner

    end type Concrete_Change_Tuner

    type, extends(Abstract_Change_Tuner), public :: Null_Change_Tuner
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: tune => Null_tune
    end type Null_Change_Tuner

contains

!implementation Abstract_Change_Tuner

    subroutine Abstract_construct(this, changed_coordinates, parameters)
        class(Abstract_Change_Tuner), intent(out) :: this
        class(Abstract_Changed_Coordinates), target, intent(in) :: changed_coordinates
        type(Concrete_Change_Tuner_Parameters), intent(in) :: parameters

        this%changed_coordinates => changed_coordinates
        call check_positive("Abstract_Change_Tuner: construct", "parameters%accumulation_period", &
            parameters%accumulation_period)
        if (num_tuning_steps < parameters%accumulation_period) then
            call warning_continue("Abstract_Change_Tuner: construct: "//&
                "num_tuning_steps < accumulation_period")
        end if
        this%accumulation_period = parameters%accumulation_period
        this%accumulated_success_ratio = 0._DP
        call check_ratio("Abstract_Change_Tuner: construct", "parameters%wanted_success_ratio", &
            parameters%wanted_success_ratio)
        call check_positive("Abstract_Change_Tuner: construct", "parameters%tolerance", &
            parameters%tolerance)
        this%success_ratio_min = parameters%wanted_success_ratio - parameters%tolerance
        this%success_ratio_max = parameters%wanted_success_ratio + parameters%tolerance
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Change_Tuner), intent(inout) :: this

        this%changed_coordinates => null()
    end subroutine Abstract_destroy

    subroutine Abstract_tune(this, tuned, i_step, success_ratio)
        class(Abstract_Change_Tuner), intent(inout) :: this
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
                call this%changed_coordinates%decrease_delta()
            else if (this%success_ratio_max < averaged_success_ratio) then
                call this%changed_coordinates%increase_delta()
            else
                tuned = .true. ! sufficient condition?
                return
            end if
            this%accumulated_success_ratio = 0._DP
        end if
    end subroutine Abstract_tune

!end implementation Abstract_Change_Tuner

!implementation Null_Change_Tuner

    subroutine Null_construct(this, changed_coordinates, parameters)
        class(Null_Change_Tuner), intent(out) :: this
        class(Abstract_Changed_Coordinates), target, intent(in) :: changed_coordinates
        type(Concrete_Change_Tuner_Parameters), intent(in) :: parameters
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Change_Tuner), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_tune(this, tuned, i_step, success_ratio)
        class(Null_Change_Tuner), intent(inout) :: this
        logical, intent(out) :: tuned
        integer, intent(in) :: i_step
        real(DP), intent(in) :: success_ratio
        tuned = .true.
    end subroutine Null_tune

!end implementation Null_Change_Tuner

end module class_change_tuner
