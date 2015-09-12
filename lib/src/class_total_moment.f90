module class_total_moment

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_dipolar_moments, only: Dipolar_Moments_Facade

implicit none

private

    type, abstract, public :: Abstract_Total_Moment
    private
        class(Dipolar_Moments_Facade), pointer :: dipolar_moments
        real(DP) :: total_moment(num_dimensions)
    contains
        procedure :: construct => Abstract_Total_Moment_construct
        procedure :: destroy => Abstract_Total_Moment_destroy
        procedure :: get => Abstract_Total_Moment_get
    end type Abstract_Total_Moment

    type, extends(Abstract_Total_Moment), public :: Concrete_Total_Moment

    end type Concrete_Total_Moment

contains

    subroutine Abstract_Total_Moment_construct(this, dipolar_moments)
        class(Abstract_Total_Moment), intent(out) :: this
        class(Dipolar_Moments_Facade), target, intent(in) :: dipolar_moments

        integer :: i_particle

        this%dipolar_moments => dipolar_moments
        this%total_moment = 0._DP
        do i_particle = 1, this%dipolar_moments%get_num()
            this%total_moment = this%total_moment + this%dipolar_moments%get(i_particle)
        end do
    end subroutine Abstract_Total_Moment_construct

    subroutine Abstract_Total_Moment_destroy(this)
        class(Abstract_Total_Moment), intent(inout) :: this

        this%dipolar_moments => null()
    end subroutine Abstract_Total_Moment_destroy

    pure function Abstract_Total_Moment_get(this) result(total_moment)
        class(Abstract_Total_Moment), intent(in) :: this
        real(DP) :: total_moment(num_dimensions)

        total_moment = this%total_moment
    end function Abstract_Total_Moment_get

end module class_total_moment
