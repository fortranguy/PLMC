module class_component_total_moment

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_3d_array
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments

implicit none

private

    type, abstract, public :: Abstract_Component_Total_Moment
    private
        class(Abstract_Component_Dipolar_Moments), pointer :: dipolar_moments
        real(DP) :: total_moment(num_dimensions)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset => Abstract_reset
        procedure :: get => Abstract_get
        procedure :: add => Abstract_add
        procedure :: remove => Abstract_remove
    end type Abstract_Component_Total_Moment

    type, extends(Abstract_Component_Total_Moment), public :: Concrete_Component_Total_Moment

    end type Concrete_Component_Total_Moment

    type, extends(Abstract_Component_Total_Moment), public :: Null_Component_Total_Moment
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: reset => Null_reset
        procedure :: get => Null_get
        procedure :: remove => Null_remove
    end type Null_Component_Total_Moment

contains

!implementation Abstract_Component_Total_Moment

    subroutine Abstract_construct(this, dipolar_moments)
        class(Abstract_Component_Total_Moment), intent(out) :: this
        class(Abstract_Component_Dipolar_Moments), target, intent(in) :: dipolar_moments

        this%dipolar_moments => dipolar_moments
        call this%reset()
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Component_Total_Moment), intent(inout) :: this

        this%dipolar_moments => null()
    end subroutine Abstract_destroy

    pure subroutine Abstract_reset(this)
        class(Abstract_Component_Total_Moment), intent(inout) :: this

        integer :: i_particle

        this%total_moment = 0._DP
        do i_particle = 1, this%dipolar_moments%get_num()
            this%total_moment = this%total_moment + this%dipolar_moments%get(i_particle)
        end do
    end subroutine Abstract_reset

    pure function Abstract_get(this) result(total_moment)
        class(Abstract_Component_Total_Moment), intent(in) :: this
        real(DP) :: total_moment(num_dimensions)

        total_moment = this%total_moment
    end function Abstract_get

    subroutine Abstract_add(this, total_moment)
        class(Abstract_Component_Total_Moment), intent(inout) :: this
        real(DP), intent(in) :: total_moment(:)

        call check_3d_array("Abstract_Component_Total_Moment: add", "total_moment", total_moment)
        this%total_moment = this%total_moment + total_moment
    end subroutine Abstract_add

    subroutine Abstract_remove(this, total_moment)
        class(Abstract_Component_Total_Moment), intent(inout) :: this
        real(DP), intent(in) :: total_moment(:)

        call check_3d_array("Abstract_Component_Total_Moment: remove", "total_moment", total_moment)
        this%total_moment = this%total_moment - total_moment
    end subroutine Abstract_remove

!end implementation Abstract_Component_Total_Moment

!implementation Null_Component_Total_Moment

    subroutine Null_construct(this, dipolar_moments)
        class(Null_Component_Total_Moment), intent(out) :: this
        class(Abstract_Component_Dipolar_Moments), target, intent(in) :: dipolar_moments
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Component_Total_Moment), intent(inout) :: this
    end subroutine Null_destroy

    pure subroutine Null_reset(this)
        class(Null_Component_Total_Moment), intent(inout) :: this
    end subroutine Null_reset

    pure function Null_get(this) result(total_moment)
        class(Null_Component_Total_Moment), intent(in) :: this
        real(DP) :: total_moment(num_dimensions)
        total_moment = 0._DP
    end function Null_get

    subroutine Null_add(this, total_moment)
        class(Null_Component_Total_Moment), intent(inout) :: this
        real(DP), intent(in) :: total_moment(:)
    end subroutine Null_add

    subroutine Null_remove(this, total_moment)
        class(Null_Component_Total_Moment), intent(inout) :: this
        real(DP), intent(in) :: total_moment(:)
    end subroutine Null_remove

!end implementation Null_Component_Total_Moment

end module class_component_total_moment
