module class_particles_total_moment

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_checks, only: check_3d_array
use class_particles_dipolar_moments, only: Abstract_Particles_Dipolar_Moments

implicit none

private

    type, abstract, public :: Abstract_Particles_Total_Moment
    private
        class(Abstract_Particles_Dipolar_Moments), pointer :: dipolar_moments
        real(DP) :: total_moment(num_dimensions)
    contains
        procedure :: construct => Abstract_Particles_Total_Moment_construct
        procedure :: destroy => Abstract_Particles_Total_Moment_destroy
        procedure :: set => Abstract_Particles_Total_Moment_set
        procedure :: reset => Abstract_Particles_Total_Moment_reset
        procedure :: get => Abstract_Particles_Total_Moment_get
    end type Abstract_Particles_Total_Moment

    type, extends(Abstract_Particles_Total_Moment), public :: Concrete_Particles_Total_Moment

    end type Concrete_Particles_Total_Moment

    type, extends(Abstract_Particles_Total_Moment), public :: Null_Particles_Total_Moment
    contains
        procedure :: construct => Null_Particles_Total_Moment_construct
        procedure :: destroy => Null_Particles_Total_Moment_destroy
        procedure :: set => Null_Particles_Total_Moment_set
        procedure :: reset => Null_Particles_Total_Moment_reset
        procedure :: get => Null_Particles_Total_Moment_get
    end type Null_Particles_Total_Moment

contains

!implementation Abstract_Particles_Total_Moment

    subroutine Abstract_Particles_Total_Moment_construct(this, dipolar_moments)
        class(Abstract_Particles_Total_Moment), intent(out) :: this
        class(Abstract_Particles_Dipolar_Moments), target, intent(in) :: dipolar_moments

        this%dipolar_moments => dipolar_moments
        call this%reset()
    end subroutine Abstract_Particles_Total_Moment_construct

    subroutine Abstract_Particles_Total_Moment_destroy(this)
        class(Abstract_Particles_Total_Moment), intent(inout) :: this

        this%dipolar_moments => null()
    end subroutine Abstract_Particles_Total_Moment_destroy

    subroutine Abstract_Particles_Total_Moment_set(this, total_moment)
        class(Abstract_Particles_Total_Moment), intent(inout) :: this
        real(DP), intent(in) :: total_moment(:)

        call check_3d_array("Abstract_Particles_Total_Moment", "total_moment", total_moment)
        this%total_moment = total_moment
    end subroutine Abstract_Particles_Total_Moment_set

    pure subroutine Abstract_Particles_Total_Moment_reset(this)
        class(Abstract_Particles_Total_Moment), intent(inout) :: this

        integer :: i_particle

        this%total_moment = 0._DP
        do i_particle = 1, this%dipolar_moments%get_num()
            this%total_moment = this%total_moment + this%dipolar_moments%get(i_particle)
        end do
    end subroutine Abstract_Particles_Total_Moment_reset

    pure function Abstract_Particles_Total_Moment_get(this) result(total_moment)
        class(Abstract_Particles_Total_Moment), intent(in) :: this
        real(DP) :: total_moment(num_dimensions)

        total_moment = this%total_moment
    end function Abstract_Particles_Total_Moment_get

!end implementation Abstract_Particles_Total_Moment

!implementation Null_Particles_Total_Moment

    subroutine Null_Particles_Total_Moment_construct(this, dipolar_moments)
        class(Null_Particles_Total_Moment), intent(out) :: this
        class(Abstract_Particles_Dipolar_Moments), target, intent(in) :: dipolar_moments
    end subroutine Null_Particles_Total_Moment_construct

    subroutine Null_Particles_Total_Moment_destroy(this)
        class(Null_Particles_Total_Moment), intent(inout) :: this
    end subroutine Null_Particles_Total_Moment_destroy

    subroutine Null_Particles_Total_Moment_set(this, total_moment)
        class(Null_Particles_Total_Moment), intent(inout) :: this
        real(DP), intent(in) :: total_moment(:)
    end subroutine Null_Particles_Total_Moment_set

    pure subroutine Null_Particles_Total_Moment_reset(this)
        class(Null_Particles_Total_Moment), intent(inout) :: this
    end subroutine Null_Particles_Total_Moment_reset

    pure function Null_Particles_Total_Moment_get(this) result(total_moment)
        class(Null_Particles_Total_Moment), intent(in) :: this
        real(DP) :: total_moment(num_dimensions)
        total_moment = 0._DP
    end function Null_Particles_Total_Moment_get

!end implementation Null_Particles_Total_Moment

end module class_particles_total_moment
