module classes_mixture_total_moment

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_array_size
use types_component_wrapper, only: Component_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Mixture_Total_Moment
    private
        type(Component_Wrapper), pointer :: components(:) => null()
        logical, allocatable :: are_dipolar(:)
        real(DP) :: total_moment(num_dimensions) = 0._DP
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset => Abstract_reset
        procedure :: is_dipolar => Abstract_is_dipolar
        procedure :: get => Abstract_get
        procedure :: add => Abstract_add
        procedure :: remove => Abstract_remove
    end type Abstract_Mixture_Total_Moment

    type, extends(Abstract_Mixture_Total_Moment), public :: Concrete_Mixture_Total_Moment

    end type Concrete_Mixture_Total_Moment

    type, extends(Abstract_Mixture_Total_Moment), public :: Null_Mixture_Total_Moment
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: reset => Null_reset
        procedure :: is_dipolar => Null_is_dipolar
        procedure :: get => Null_get
        procedure :: add => Null_add
        procedure :: remove => Null_remove
    end type Null_Mixture_Total_Moment

contains

!implementation Abstract_Mixture_Total_Moment

    subroutine Abstract_construct(this, components, are_dipolar)
        class(Abstract_Mixture_Total_Moment), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)

        this%components => components
        allocate(this%are_dipolar(size(are_dipolar)))
        this%are_dipolar = are_dipolar
        call this%reset()
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Mixture_Total_Moment), intent(inout) :: this

        if (allocated(this%are_dipolar)) deallocate(this%are_dipolar)
        this%components => null()
    end subroutine Abstract_destroy

    pure subroutine Abstract_reset(this)
        class(Abstract_Mixture_Total_Moment), intent(inout) :: this

        integer :: i_component, i_particle

        this%total_moment = 0._DP
        do i_component = 1, size(this%components)
            do i_particle = 1, this%components(i_component)%dipole_moments%get_num()
                this%total_moment = this%total_moment + this%components(i_component)%&
                    dipole_moments%get(i_particle)
            end do
        end do
    end subroutine Abstract_reset

    pure logical function Abstract_is_dipolar(this, i_component) result(is_dipolar)
        class(Abstract_Mixture_Total_Moment), intent(in) :: this
        integer, intent(in) :: i_component

        is_dipolar = this%are_dipolar(i_component)
    end function Abstract_is_dipolar

    pure function Abstract_get(this) result(total_moment)
        class(Abstract_Mixture_Total_Moment), intent(in) :: this
        real(DP) :: total_moment(num_dimensions)

        total_moment = this%total_moment
    end function Abstract_get

    subroutine Abstract_add(this, i_component, dipole_moment)
        class(Abstract_Mixture_Total_Moment), intent(inout) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: dipole_moment(:)

        if (.not.this%are_dipolar(i_component)) return

        call check_array_size("Abstract_Mixture_Total_Moment: add", "dipole_moment", &
            dipole_moment, num_dimensions)
        this%total_moment = this%total_moment + dipole_moment
    end subroutine Abstract_add

    subroutine Abstract_remove(this, i_component, dipole_moment)
        class(Abstract_Mixture_Total_Moment), intent(inout) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: dipole_moment(:)

        if (.not.this%are_dipolar(i_component)) return

        call check_array_size("Abstract_Mixture_Total_Moment: remove", "dipole_moment", &
            dipole_moment, num_dimensions)
        this%total_moment = this%total_moment - dipole_moment
    end subroutine Abstract_remove

!end implementation Abstract_Mixture_Total_Moment

!implementation Null_Mixture_Total_Moment

    subroutine Null_construct(this, components, are_dipolar)
        class(Null_Mixture_Total_Moment), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Mixture_Total_Moment), intent(inout) :: this
    end subroutine Null_destroy

    pure function Null_get(this) result(total_moment)
        class(Null_Mixture_Total_Moment), intent(in) :: this
        real(DP) :: total_moment(num_dimensions)
        total_moment = 0._DP
    end function Null_get

    pure subroutine Null_reset(this)
        class(Null_Mixture_Total_Moment), intent(inout) :: this
    end subroutine Null_reset

    pure logical function Null_is_dipolar(this, i_component) result(is_dipolar)
        class(Null_Mixture_Total_Moment), intent(in) :: this
        integer, intent(in) :: i_component
        is_dipolar = .false.
    end function Null_is_dipolar

    subroutine Null_add(this, i_component, dipole_moment)
        class(Null_Mixture_Total_Moment), intent(inout) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: dipole_moment(:)
    end subroutine Null_add

    subroutine Null_remove(this, i_component, dipole_moment)
        class(Null_Mixture_Total_Moment), intent(inout) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: dipole_moment(:)
    end subroutine Null_remove

!end implementation Null_Mixture_Total_Moment

end module classes_mixture_total_moment
