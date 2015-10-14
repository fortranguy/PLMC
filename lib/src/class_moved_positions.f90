module class_moved_positions

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_errors, only: warning_continue
use procedures_checks, only: check_3d_array, check_positive, check_increase_factor
use class_periodic_box, only: Abstract_Periodic_Box
use class_component_positions, only: Abstract_Component_Positions
use module_change_tuning, only: Concrete_Tuning_Parameters, &
    set_increase_factor
use class_changed_coordinates, only: Abstract_Changed_Coordinates

implicit none

private

    type, extends(Abstract_Changed_Coordinates), abstract, public :: Abstract_Moved_Positions
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Component_Positions), pointer :: positions => null()
        real(DP) :: delta(num_dimensions)
        type(Concrete_Tuning_Parameters) :: tuning_parameters
        real(DP) :: current_increase_factor
        logical :: max_factor_reached
    contains
        procedure :: construct => Abstract_Moved_Positions_construct
        procedure :: destroy => Abstract_Moved_Positions_destroy
        procedure :: increase_delta => Abstract_Moved_Positions_increase_delta
        procedure :: decrease_delta => Abstract_Moved_Positions_decrease_delta
        procedure :: get => Abstract_Moved_Positions_get
    end type Abstract_Moved_Positions

    type, extends(Abstract_Moved_Positions), public :: Concrete_Moved_Positions

    end type Concrete_Moved_Positions

    type, extends(Abstract_Moved_Positions), public :: Null_Moved_Positions
    contains
        procedure :: construct => Null_Moved_Positions_construct
        procedure :: destroy => Null_Moved_Positions_destroy
        procedure :: increase_delta => Null_Moved_Positions_increase_delta
        procedure :: decrease_delta => Null_Moved_Positions_decrease_delta
        procedure :: get => Null_Moved_Positions_get
    end type Null_Moved_Positions

contains

!implementation Abstract_Moved_Positions

    subroutine Abstract_Moved_Positions_construct(this, periodic_box, positions, delta, &
        tuning_parameters)
        class(Abstract_Moved_Positions), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Positions), target, intent(in) :: positions
        real(DP), intent(in) :: delta(:)
        type(Concrete_Tuning_Parameters), intent(in) :: tuning_parameters

        this%periodic_box => periodic_box
        this%positions => positions
        call check_3d_array("Abstract_Moved_Positions", "delta", delta)
        call check_positive("Abstract_Moved_Positions", "delta", delta)
        this%delta = delta
        call check_increase_factor("Abstract_Moved_Positions", "increase_factor", &
            tuning_parameters%increase_factor)
        this%tuning_parameters%increase_factor = tuning_parameters%increase_factor
        this%current_increase_factor = this%tuning_parameters%increase_factor
        call check_increase_factor("Abstract_Moved_Positions", "increase_factor_max", &
            tuning_parameters%increase_factor_max)
        this%tuning_parameters%increase_factor_max = tuning_parameters%increase_factor_max
    end subroutine Abstract_Moved_Positions_construct

    subroutine Abstract_Moved_Positions_destroy(this)
        class(Abstract_Moved_Positions), intent(inout) :: this

        this%positions => null()
        this%periodic_box => null()
    end subroutine Abstract_Moved_Positions_destroy

    subroutine Abstract_Moved_Positions_increase_delta(this)
        class(Abstract_Moved_Positions), intent(inout) :: this

        if (this%max_factor_reached) return
        call set_increase_factor("Abstract_Moved_Positions", this%current_increase_factor, &
            this%tuning_parameters, this%max_factor_reached)
        this%delta = this%current_increase_factor * this%delta
    end subroutine Abstract_Moved_Positions_increase_delta

    subroutine Abstract_Moved_Positions_decrease_delta(this)
        class(Abstract_Moved_Positions), intent(inout) :: this

        this%delta = this%delta / this%current_increase_factor
    end subroutine Abstract_Moved_Positions_decrease_delta

    function Abstract_Moved_Positions_get(this, i_particle) result(moved_position)
        real(DP) :: moved_position(num_dimensions)
        class(Abstract_Moved_Positions), intent(in) :: this
        integer, intent(in) :: i_particle

        real(DP) :: rand(num_dimensions)

        call random_number(rand)
        moved_position = this%positions%get(i_particle) + (rand - 0.5_DP) * this%delta
        moved_position = this%periodic_box%folded(moved_position)
    end function Abstract_Moved_Positions_get

!end implementation Abstract_Moved_Positions

!implementation Null_Moved_Positions

    subroutine Null_Moved_Positions_construct(this, periodic_box, positions, delta, &
        tuning_parameters)
        class(Null_Moved_Positions), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Positions), target, intent(in) :: positions
        real(DP), intent(in) :: delta(:)
        type(Concrete_Tuning_Parameters), intent(in) :: tuning_parameters
    end subroutine Null_Moved_Positions_construct

    subroutine Null_Moved_Positions_destroy(this)
        class(Null_Moved_Positions), intent(inout) :: this
    end subroutine Null_Moved_Positions_destroy

    subroutine Null_Moved_Positions_increase_delta(this)
        class(Null_Moved_Positions), intent(inout) :: this
    end subroutine Null_Moved_Positions_increase_delta

    subroutine Null_Moved_Positions_decrease_delta(this)
        class(Null_Moved_Positions), intent(inout) :: this
    end subroutine Null_Moved_Positions_decrease_delta

    function Null_Moved_Positions_get(this, i_particle) result(moved_position)
        class(Null_Moved_Positions), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: moved_position(num_dimensions)
        moved_position = 0._DP
    end function Null_Moved_Positions_get

!end implementation Null_Moved_Positions

end module class_moved_positions
