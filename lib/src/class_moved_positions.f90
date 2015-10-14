module class_moved_positions

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_errors, only: warning_continue
use procedures_checks, only: check_3d_array, check_positive, check_increase_factor
use class_periodic_box, only: Abstract_Periodic_Box
use class_component_coordinates, only: Abstract_Component_Coordinates
use module_change_tuning, only: Concrete_Tuning_Parameters, &
    set_increase_factor
use class_changed_coordinates, only: Abstract_Changed_Coordinates

implicit none

private

    type, extends(Abstract_Changed_Coordinates), public :: Concrete_Moved_Positions
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Component_Coordinates), pointer :: positions => null()
        real(DP) :: delta(num_dimensions)
        type(Concrete_Tuning_Parameters) :: tuning_parameters
        real(DP) :: current_increase_factor
        logical :: max_factor_reached
    contains
        procedure :: construct => Concrete_Moved_Positions_construct
        procedure :: destroy => Concrete_Moved_Positions_destroy
        procedure :: increase_delta => Concrete_Moved_Positions_increase_delta
        procedure :: decrease_delta => Concrete_Moved_Positions_decrease_delta
        procedure :: get => Concrete_Moved_Positions_get
    end type Concrete_Moved_Positions

contains

!implementation Concrete_Moved_Positions

    subroutine Concrete_Moved_Positions_construct(this, periodic_box, positions, delta, &
        tuning_parameters)
        class(Concrete_Moved_Positions), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: positions
        real(DP), intent(in) :: delta(:)
        type(Concrete_Tuning_Parameters), intent(in) :: tuning_parameters

        this%periodic_box => periodic_box
        this%positions => positions
        call check_3d_array("Concrete_Moved_Positions", "delta", delta)
        call check_positive("Concrete_Moved_Positions", "delta", delta)
        this%delta = delta
        call check_increase_factor("Concrete_Moved_Positions", "increase_factor", &
            tuning_parameters%increase_factor)
        this%tuning_parameters%increase_factor = tuning_parameters%increase_factor
        this%current_increase_factor = this%tuning_parameters%increase_factor
        call check_increase_factor("Concrete_Moved_Positions", "increase_factor_max", &
            tuning_parameters%increase_factor_max)
        this%tuning_parameters%increase_factor_max = tuning_parameters%increase_factor_max
    end subroutine Concrete_Moved_Positions_construct

    subroutine Concrete_Moved_Positions_destroy(this)
        class(Concrete_Moved_Positions), intent(inout) :: this

        this%positions => null()
        this%periodic_box => null()
    end subroutine Concrete_Moved_Positions_destroy

    subroutine Concrete_Moved_Positions_increase_delta(this)
        class(Concrete_Moved_Positions), intent(inout) :: this

        if (this%max_factor_reached) return
        call set_increase_factor("Concrete_Moved_Positions", this%current_increase_factor, &
            this%tuning_parameters, this%max_factor_reached)
        this%delta = this%current_increase_factor * this%delta
    end subroutine Concrete_Moved_Positions_increase_delta

    subroutine Concrete_Moved_Positions_decrease_delta(this)
        class(Concrete_Moved_Positions), intent(inout) :: this

        this%delta = this%delta / this%current_increase_factor
    end subroutine Concrete_Moved_Positions_decrease_delta

    function Concrete_Moved_Positions_get(this, i_particle) result(moved_position)
        real(DP) :: moved_position(num_dimensions)
        class(Concrete_Moved_Positions), intent(in) :: this
        integer, intent(in) :: i_particle

        real(DP) :: rand(num_dimensions)

        call random_number(rand)
        moved_position = this%positions%get(i_particle) + (rand - 0.5_DP) * this%delta
        moved_position = this%periodic_box%folded(moved_position)
    end function Concrete_Moved_Positions_get

!end implementation Concrete_Moved_Positions

end module class_moved_positions
