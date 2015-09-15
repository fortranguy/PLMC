module class_total_moment

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_checks, only: check_3d_array
use class_dipolar_moments, only: Dipolar_Moments_Facade

implicit none

private

    type, public :: Total_Moment_Facade
    private
        class(Dipolar_Moments_Facade), pointer :: dipolar_moments
        real(DP) :: total_moment(num_dimensions)
    contains
        procedure :: construct => Total_Moment_Facade_construct
        procedure :: destroy => Total_Moment_Facade_destroy
        procedure :: set => Total_Moment_Facade_set
        procedure :: reset => Total_Moment_Facade_reset
        procedure :: get => Total_Moment_Facade_get
    end type Total_Moment_Facade

contains

    subroutine Total_Moment_Facade_construct(this, dipolar_moments)
        class(Total_Moment_Facade), intent(out) :: this
        class(Dipolar_Moments_Facade), target, intent(in) :: dipolar_moments

        this%dipolar_moments => dipolar_moments
        call this%reset()
    end subroutine Total_Moment_Facade_construct

    subroutine Total_Moment_Facade_destroy(this)
        class(Total_Moment_Facade), intent(inout) :: this

        this%dipolar_moments => null()
    end subroutine Total_Moment_Facade_destroy

    subroutine Total_Moment_Facade_set(this, total_moment)
        class(Total_Moment_Facade), intent(inout) :: this
        real(DP), intent(in) :: total_moment(:)

        call check_3d_array("Total_Moment_Facade", "total_moment", total_moment)
        this%total_moment = total_moment
    end subroutine Total_Moment_Facade_set

    pure subroutine Total_Moment_Facade_reset(this)
        class(Total_Moment_Facade), intent(inout) :: this

        integer :: i_particle

        this%total_moment = 0._DP
        do i_particle = 1, this%dipolar_moments%get_num()
            this%total_moment = this%total_moment + this%dipolar_moments%get(i_particle)
        end do
    end subroutine Total_Moment_Facade_reset

    pure function Total_Moment_Facade_get(this) result(total_moment)
        class(Total_Moment_Facade), intent(in) :: this
        real(DP) :: total_moment(num_dimensions)

        total_moment = this%total_moment
    end function Total_Moment_Facade_get

end module class_total_moment
