module class_moved_positions

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_errors, only: warning_continue
use procedures_checks, only: check_3d_array, check_positive, check_increase_factor
use class_particles_positions, only: Abstract_Particles_Positions
use module_adaptation, only: Concrete_Adaptation_Parameters, &
    set_increase_factor

implicit none

private

    type, abstract, public :: Abstract_Moved_Positions
    private
        class(Abstract_Particles_Positions), pointer :: positions
        real(DP) :: delta(num_dimensions)
        type(Concrete_Adaptation_Parameters) :: adaptation_parameters
        real(DP) :: current_increase_factor
        logical :: max_factor_reached
    contains
        procedure :: construct => Abstract_Moved_Positions_construct
        procedure :: destroy => Abstract_Moved_Positions_destroy
        procedure :: increase => Abstract_Moved_Positions_increase
        procedure :: decrease => Abstract_Moved_Positions_decrease
        procedure :: get => Abstract_Moved_Positions_get
    end type Abstract_Moved_Positions

    type, extends(Abstract_Moved_Positions), public :: Concrete_Moved_Positions

    end type Concrete_Moved_Positions

contains

    subroutine Abstract_Moved_Positions_construct(this, positions, delta, adaptation_parameters)
        class(Abstract_Moved_Positions), intent(out) :: this
        class(Abstract_Particles_Positions), target, intent(in) :: positions
        real(DP), intent(in) :: delta(:)
        type(Concrete_Adaptation_Parameters) :: adaptation_parameters

        this%positions => positions
        call check_3d_array("Abstract_Moved_Positions", "delta", delta)
        call check_positive("Abstract_Moved_Positions", "delta", delta)
        this%delta = delta
        call check_increase_factor("Abstract_Moved_Positions", "increase_factor", &
            adaptation_parameters%increase_factor)
        this%adaptation_parameters%increase_factor = adaptation_parameters%increase_factor
        this%current_increase_factor = this%adaptation_parameters%increase_factor
        call check_increase_factor("Abstract_Moved_Positions", "increase_factor_max", &
            adaptation_parameters%increase_factor_max)
        this%adaptation_parameters%increase_factor_max = adaptation_parameters%increase_factor_max
    end subroutine Abstract_Moved_Positions_construct

    subroutine Abstract_Moved_Positions_destroy(this)
        class(Abstract_Moved_Positions), intent(inout) :: this

        this%positions => null()
    end subroutine Abstract_Moved_Positions_destroy

    subroutine Abstract_Moved_Positions_increase(this)
        class(Abstract_Moved_Positions), intent(inout) :: this

        if (this%max_factor_reached) return
        call set_increase_factor("Abstract_Moved_Positions", this%current_increase_factor, &
            this%adaptation_parameters, this%max_factor_reached)
        this%delta = this%current_increase_factor * this%delta
    end subroutine Abstract_Moved_Positions_increase

    subroutine Abstract_Moved_Positions_decrease(this)
        class(Abstract_Moved_Positions), intent(inout) :: this

        this%delta = this%delta / this%current_increase_factor
    end subroutine Abstract_Moved_Positions_decrease

    function Abstract_Moved_Positions_get(this, i_particle) result(moved_position)
        class(Abstract_Moved_Positions), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: moved_position(num_dimensions)

        real(DP) :: rand(num_dimensions)

        call random_number(rand)
        moved_position = this%positions%get(i_particle) + (rand - 0.5_DP) * this%delta
    end function Abstract_Moved_Positions_get

end module class_moved_positions
