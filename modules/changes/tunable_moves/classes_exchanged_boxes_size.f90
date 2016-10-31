module classes_exchanged_boxes_size

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive, check_increase_factor
use classes_tunable_move, only: Abstract_Tunable_Move
use module_move_tuning, only: Concrete_Move_Tuning_Parameters, set_increase_factor

implicit none

private

    type, extends(Abstract_Tunable_Move), abstract, public :: Abstract_Exchanged_Boxes_Size
    private
        real(DP) :: frequency_ratio = 0._DP
        real(DP) :: delta = 0._DP
        type(Concrete_Move_Tuning_Parameters) :: tuning_parameters
        real(DP) :: current_increase_factor = 1._DP
        logical :: max_factor_reached = .false.
    contains
        procedure :: set => Abstract_set
        procedure :: get_frequency_ratio => Abstract_get_frequency_ratio
        procedure(Abstract_get_ratios), deferred :: get_ratios
        procedure :: increase_delta => Abstract_increase_delta
        procedure :: decrease_delta => Abstract_decrease_delta
    end type Abstract_Exchanged_Boxes_Size

    type, extends(Abstract_Exchanged_Boxes_Size), public :: XYZ_Exchanged_Boxes_Size
    contains
        procedure :: get_ratios => XYZ_get_ratios
    end type XYZ_Exchanged_Boxes_Size

    type, extends(Abstract_Exchanged_Boxes_Size), public :: XY_Exchanged_Boxes_Size
    contains
        procedure :: get_ratios => XY_get_ratios
    end type XY_Exchanged_Boxes_Size

    type, extends(Abstract_Exchanged_Boxes_Size), public :: Null_Exchanged_Boxes_Size
    contains
        procedure :: set => Null_set
        procedure :: get_frequency_ratio => Null_get_frequency_ratio
        procedure :: get_ratios => Null_get_ratios
        procedure :: increase_delta => Null_increase_delta
        procedure :: decrease_delta => Null_decrease_delta
    end type Null_Exchanged_Boxes_Size

    abstract interface

        !> Boxes size exchange such that:
        !> \[
        !>      \ln \left( \frac{V_{\boldsymbol{I}}^\prime}{V_{\boldsymbol{J}}^\prime} \right) =
        !>          \ln \left( \frac{V_{\boldsymbol{I}}}{V_{\boldsymbol{J}}} \right) +
        !>          \mathrm{rand}[-1/2, +1/2] \delta
        !> \]
        function Abstract_get_ratios(this, hetero_boxes_size_ratio) result(ratios)
        import :: DP, num_dimensions, Abstract_Exchanged_Boxes_Size
            real(DP) :: ratios(num_dimensions, 2)
            class(Abstract_Exchanged_Boxes_Size), intent(in) :: this
            real(DP), intent(in) :: hetero_boxes_size_ratio(:)
                !! \( \vec{L}_{\boldsymbol{I}} / \vec{L}_{\boldsymbol{J}} \)
        end function Abstract_get_ratios

    end interface

contains

!implementation Abstract_Exchanged_Boxes_Size

    subroutine Abstract_set(this, frequency_ratio, initial_delta, tuning_parameters)
        class(Abstract_Exchanged_Boxes_Size), intent(inout) :: this
        real(DP), intent(in) :: frequency_ratio
        real(DP), intent(in) :: initial_delta
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters

        call check_positive("Abstract_Exchanged_Boxes_Size: construct", "frequency_ratio", &
            frequency_ratio)
        this%frequency_ratio = frequency_ratio
        call check_positive("Abstract_Exchanged_Boxes_Size: construct", "initial_delta", &
            initial_delta)
        this%delta = initial_delta

        call check_increase_factor("Abstract_Exchanged_Boxes_Size: construct", "increase_factor", &
            tuning_parameters%increase_factor)
        this%tuning_parameters%increase_factor = tuning_parameters%increase_factor
        this%current_increase_factor = this%tuning_parameters%increase_factor
        call check_increase_factor("Abstract_Exchanged_Boxes_Size: construct", &
            "increase_factor_max", tuning_parameters%increase_factor_max)
        this%tuning_parameters%increase_factor_max = tuning_parameters%increase_factor_max
    end subroutine Abstract_set

    pure real(DP) function Abstract_get_frequency_ratio(this) result(frequency_ratio)
        class(Abstract_Exchanged_Boxes_Size), intent(in) :: this

        frequency_ratio = this%frequency_ratio
    end function Abstract_get_frequency_ratio

    subroutine Abstract_increase_delta(this)
        class(Abstract_Exchanged_Boxes_Size), intent(inout) :: this

        if (this%max_factor_reached) return
        call set_increase_factor("Abstract_Exchanged_Boxes_Size", this%current_increase_factor, &
            this%tuning_parameters, this%max_factor_reached)
        this%delta = this%current_increase_factor * this%delta
    end subroutine Abstract_increase_delta

    subroutine Abstract_decrease_delta(this)
        class(Abstract_Exchanged_Boxes_Size), intent(inout) :: this

        this%delta = this%delta / this%current_increase_factor
    end subroutine Abstract_decrease_delta

!end implementation Abstract_Exchanged_Boxes_Size

!implementation XYZ_Exchanged_Boxes_Size

    !> \[
    !>      \frac{V_{\boldsymbol{I}}^\prime}{V_{\boldsymbol{I}}} =
    !>          \frac{\left( 1 + \frac{V_{\boldsymbol{I}}}{V_{\boldsymbol{J}}} \right)
    !>              e^{\mathrm{rand}[-1/2, +1/2] \delta}}
    !>              {1 + \frac{V_{\boldsymbol{I}}}
    !>                  {V_{\boldsymbol{J}}}e^{\mathrm{rand}[-1/2, +1/2] \delta}}
    !> \]
    !> \[
    !>      \frac{V_{\boldsymbol{J}}^\prime}{V_{\boldsymbol{J}}} =
    !>          \frac{V_{\boldsymbol{I}}^\prime}
    !>              {V_{\boldsymbol{I}}}e^{-\mathrm{rand}[-1/2, +1/2] \delta}
    !> \]
    function XYZ_get_ratios(this, hetero_boxes_size_ratio) result(ratios)
        real(DP) :: ratios(num_dimensions, 2)
        class(XYZ_Exchanged_Boxes_Size), intent(in) :: this
        real(DP), intent(in) :: hetero_boxes_size_ratio(:)

        real(DP) :: exp_rand_delta, rand

        call random_number(rand)
        exp_rand_delta = exp((rand - 0.5_DP) * this%delta)

        ratios(:, 1) = (1._DP + product(hetero_boxes_size_ratio)) * exp_rand_delta / &
            (1._DP + product(hetero_boxes_size_ratio) * exp_rand_delta)
        ratios(:, 2) = ratios(:, 1) / exp_rand_delta
        ratios = ratios ** (1._DP/3._DP)
    end function XYZ_get_ratios

!end implementation XYZ_Exchanged_Boxes_Size

!implementation XY_Exchanged_Boxes_Size

    !> \[
    !>      \frac{S_{\boldsymbol{I}}^\prime}{S_{\boldsymbol{I}}} =
    !>          \frac{\left( 1 + \frac{S_{\boldsymbol{I}}}{S_{\boldsymbol{J}}} \right)
    !>              e^{\mathrm{rand}[-1/2, +1/2] \delta}}
    !>              {1 + \frac{S_{\boldsymbol{I}}}{S_{\boldsymbol{J}}}
    !>                  e^{\mathrm{rand}[-1/2, +1/2] \delta}}
    !> \]
    !> \[
    !>      \frac{S_{\boldsymbol{J}}^\prime}{S_{\boldsymbol{J}}} =
    !>          \frac{S_{\boldsymbol{I}}^\prime}{S_{\boldsymbol{I}}}
    !>              e^{-\mathrm{rand}[-1/2, +1/2] \delta}
    !> \]
    function XY_get_ratios(this, hetero_boxes_size_ratio) result(ratios)
        real(DP) :: ratios(num_dimensions, 2)
        class(XY_Exchanged_Boxes_Size), intent(in) :: this
        real(DP), intent(in) :: hetero_boxes_size_ratio(:)

        real(DP) :: exp_rand_delta, rand

        call random_number(rand)
        exp_rand_delta = exp((rand - 0.5_DP) * this%delta)

        ratios(1:2, 1) = (1._DP + product(hetero_boxes_size_ratio(1:2))) * exp_rand_delta / &
            (1._DP + product(hetero_boxes_size_ratio(1:2)) * exp_rand_delta)
        ratios(1:2, 2) = ratios(1:2, 1) / exp_rand_delta

        ratios(1:2, :) = ratios(1:2, :) ** (1._DP/2._DP)
        ratios(3, :) = 1._DP
    end function XY_get_ratios

!end implementation XY_Exchanged_Boxes_Size

!implementation Null_Exchanged_Boxes_Size

    subroutine Null_set(this, frequency_ratio, initial_delta, tuning_parameters)
        class(Null_Exchanged_Boxes_Size), intent(inout) :: this
        real(DP), intent(in) :: frequency_ratio
        real(DP), intent(in) :: initial_delta
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
    end subroutine Null_set

    pure real(DP) function Null_get_frequency_ratio(this) result(frequency_ratio)
        class(Null_Exchanged_Boxes_Size), intent(in) :: this
        frequency_ratio = 0._DP
    end function Null_get_frequency_ratio

    function Null_get_ratios(this, hetero_boxes_size_ratio) result(ratios)
        real(DP) :: ratios(num_dimensions, 2)
        class(Null_Exchanged_Boxes_Size), intent(in) :: this
        real(DP), intent(in) :: hetero_boxes_size_ratio(:)
        ratios = 1._DP
    end function Null_get_ratios

    subroutine Null_increase_delta(this)
        class(Null_Exchanged_Boxes_Size), intent(inout) :: this
    end subroutine Null_increase_delta

    subroutine Null_decrease_delta(this)
        class(Null_Exchanged_Boxes_Size), intent(inout) :: this
    end subroutine Null_decrease_delta

!end implementation Null_Exchanged_Boxes_Size

end module classes_exchanged_boxes_size
