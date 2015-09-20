module class_particles_total_moment

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_checks, only: check_3d_array
use class_particles_dipolar_moments, only: Particles_Dipolar_Moments_Facade

implicit none

private

    type, public :: Particles_Total_Moment_Facade
    private
        class(Particles_Dipolar_Moments_Facade), pointer :: dipolar_moments
        real(DP) :: total_moment(num_dimensions)
    contains
        procedure :: construct => Particles_Total_Moment_Facade_construct
        procedure :: destroy => Particles_Total_Moment_Facade_destroy
        procedure :: set => Particles_Total_Moment_Facade_set
        procedure :: reset => Particles_Total_Moment_Facade_reset
        procedure :: get => Particles_Total_Moment_Facade_get
    end type Particles_Total_Moment_Facade

contains

    subroutine Particles_Total_Moment_Facade_construct(this, dipolar_moments)
        class(Particles_Total_Moment_Facade), intent(out) :: this
        class(Particles_Dipolar_Moments_Facade), target, intent(in) :: dipolar_moments

        this%dipolar_moments => dipolar_moments
        call this%reset()
    end subroutine Particles_Total_Moment_Facade_construct

    subroutine Particles_Total_Moment_Facade_destroy(this)
        class(Particles_Total_Moment_Facade), intent(inout) :: this

        this%dipolar_moments => null()
    end subroutine Particles_Total_Moment_Facade_destroy

    subroutine Particles_Total_Moment_Facade_set(this, total_moment)
        class(Particles_Total_Moment_Facade), intent(inout) :: this
        real(DP), intent(in) :: total_moment(:)

        call check_3d_array("Particles_Total_Moment_Facade", "total_moment", total_moment)
        this%total_moment = total_moment
    end subroutine Particles_Total_Moment_Facade_set

    pure subroutine Particles_Total_Moment_Facade_reset(this)
        class(Particles_Total_Moment_Facade), intent(inout) :: this

        integer :: i_particle

        this%total_moment = 0._DP
        do i_particle = 1, this%dipolar_moments%get_num()
            this%total_moment = this%total_moment + this%dipolar_moments%get(i_particle)
        end do
    end subroutine Particles_Total_Moment_Facade_reset

    pure function Particles_Total_Moment_Facade_get(this) result(total_moment)
        class(Particles_Total_Moment_Facade), intent(in) :: this
        real(DP) :: total_moment(num_dimensions)

        total_moment = this%total_moment
    end function Particles_Total_Moment_Facade_get

end module class_particles_total_moment
