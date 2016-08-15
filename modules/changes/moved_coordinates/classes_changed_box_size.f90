module classes_changed_box_size

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive, check_increase_factor
use classes_moved_coordinates, only: Abstract_Moved_Coordinates
use module_move_tuning, only: Concrete_Move_Tuning_Parameters, set_increase_factor

implicit none

private

    type, extends(Abstract_Moved_Coordinates), abstract, public :: Abstract_Changed_Box_Size
    private
        real(DP) :: delta = 0._DP
        type(Concrete_Move_Tuning_Parameters) :: tuning_parameters
        real(DP) :: current_increase_factor = 1._DP
        logical :: max_factor_reached = .false.
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: increase_delta => Abstract_increase_delta
        procedure :: decrease_delta => Abstract_decrease_delta
        procedure(Abstract_get_ratio), deferred :: get_ratio
    end type Abstract_Changed_Box_Size

    abstract interface

        !> Box size change such that:
        !> \[ \frac{V^\prime}{V} = e^{\mathrm{rand}[-1/2, +1/2] \delta} \]
        function Abstract_get_ratio(this) result(ratio)
        import :: DP, num_dimensions, Abstract_Changed_Box_Size
            class(Abstract_Changed_Box_Size), intent(in) :: this
            real(DP) :: ratio(num_dimensions)
        end function Abstract_get_ratio

    end interface

    type, extends(Abstract_Changed_Box_Size), public :: XYZ_Changed_Box_Size
    contains
        procedure :: get_ratio => XYZ_get_ratio
    end type XYZ_Changed_Box_Size

    type, extends(Abstract_Changed_Box_Size), public :: XY_Changed_Box_Size
    contains
        procedure :: get_ratio => XY_get_ratio
    end type

    type, extends(Abstract_Changed_Box_Size), public :: Null_Changed_Box_Size
    contains
        procedure :: construct => Null_construct
        procedure :: increase_delta => Null_increase_delta
        procedure :: decrease_delta => Null_decrease_delta
        procedure :: get_ratio => Null_get_ratio
    end type Null_Changed_Box_Size

contains

!implementation Abstract_Changed_Box_Size

    subroutine Abstract_construct(this, initial_delta, tuning_parameters)
        class(Abstract_Changed_Box_Size), intent(out) :: this
        real(DP), intent(in) :: initial_delta
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters

        call check_positive("Abstract_Changed_Box_Size: construct", "initial_delta", initial_delta)
        this%delta = initial_delta
        call check_increase_factor("Abstract_Changed_Box_Size: construct", "increase_factor", &
            tuning_parameters%increase_factor)
        this%tuning_parameters%increase_factor = tuning_parameters%increase_factor
        this%current_increase_factor = this%tuning_parameters%increase_factor
        call check_increase_factor("Abstract_Changed_Box_Size: construct", &
            "increase_factor_max", tuning_parameters%increase_factor_max)
        this%tuning_parameters%increase_factor_max = tuning_parameters%increase_factor_max
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Changed_Box_Size), intent(inout) :: this
    end subroutine Abstract_destroy

    subroutine Abstract_increase_delta(this)
        class(Abstract_Changed_Box_Size), intent(inout) :: this

        if (this%max_factor_reached) return
        call set_increase_factor("Abstract_Changed_Box_Size", this%current_increase_factor, &
            this%tuning_parameters, this%max_factor_reached)
        this%delta = this%current_increase_factor * this%delta
    end subroutine Abstract_increase_delta

    subroutine Abstract_decrease_delta(this)
        class(Abstract_Changed_Box_Size), intent(inout) :: this

        this%delta = this%delta / this%current_increase_factor
    end subroutine Abstract_decrease_delta

!end implementation Abstract_Changed_Box_Size

!implementation XYZ_Changed_Box_Size

    !> \[ \frac{L_{1:3}^\prime}{L_{1:3}} = e^{\mathrm{rand}[-1/2, +1/2] \delta / 3} \]
    function XYZ_get_ratio(this) result(ratio)
        real(DP) :: ratio(num_dimensions)
        class(XYZ_Changed_Box_Size), intent(in) :: this

        real(DP) :: rand

        call random_number(rand)
        ratio = exp((rand - 0.5_DP) * this%delta / 3._DP)
    end function XYZ_get_ratio

!end implementation XYZ_Changed_Box_Size

!implementation XY_Changed_Box_Size

    !> \[
    !>      \left( \frac{L_{1:2}^\prime}{L_{1:2}}, \frac{L_3^\prime}{L_3} \right) =
    !>          \left( e^{\mathrm{rand}[-1/2, +1/2] \delta / 2}, 1 \right)
    !> \]
    function XY_get_ratio(this) result(ratio)
        real(DP) :: ratio(num_dimensions)
        class(XY_Changed_Box_Size), intent(in) :: this

        real(DP) :: rand

        call random_number(rand)
        ratio(1:2) = exp((rand - 0.5_DP) * this%delta / 2._DP)
        ratio(3) = 1._DP
    end function XY_get_ratio

!end implementation XY_Changed_Box_Size

!implementation Null_Changed_Box_Size

    subroutine Null_construct(this, initial_delta, tuning_parameters)
        class(Null_Changed_Box_Size), intent(out) :: this
        real(DP), intent(in) :: initial_delta
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
    end subroutine Null_construct

    subroutine Null_increase_delta(this)
        class(Null_Changed_Box_Size), intent(inout) :: this
    end subroutine Null_increase_delta

    subroutine Null_decrease_delta(this)
        class(Null_Changed_Box_Size), intent(inout) :: this
    end subroutine Null_decrease_delta

    function Null_get_ratio(this) result(ratio)
        real(DP) :: ratio(num_dimensions)
        class(Null_Changed_Box_Size), intent(in) :: this
        ratio = 1._DP
    end function Null_get_ratio

!end implementation Null_Changed_Box_Size


end module classes_changed_box_size
