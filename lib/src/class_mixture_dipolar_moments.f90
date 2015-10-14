module class_mixture_dipolar_moments

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use class_coordinates, only: Abstract_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments

implicit none

private

    type, extends(Abstract_Coordinates), abstract, public :: Abstract_Mixture_Dipolar_Moments
    private
        class(Abstract_Component_Dipolar_Moments), pointer :: component_1_moments => null(), &
            component_2_moments => null()
    contains
        procedure :: construct => Abstract_Mixture_Dipolar_Moments_construct
        procedure :: destroy => Abstract_Mixture_Dipolar_Moments_destroy
        procedure :: get_num => Abstract_Mixture_Dipolar_Moments_get_num
        procedure :: get => Abstract_Mixture_Dipolar_Moments_get
    end type Abstract_Mixture_Dipolar_Moments

contains

    subroutine Abstract_Mixture_Dipolar_Moments_construct(this, component_1_moments, &
        component_2_moments)
        class(Abstract_Mixture_Dipolar_Moments), intent(out) :: this
        class(Abstract_Component_Dipolar_Moments), target, intent(in) :: component_1_moments, &
            component_2_moments

        this%component_1_moments => component_1_moments
        this%component_2_moments => component_2_moments
    end subroutine Abstract_Mixture_Dipolar_Moments_construct

    subroutine Abstract_Mixture_Dipolar_Moments_destroy(this)
        class(Abstract_Mixture_Dipolar_Moments), intent(inout) :: this

        this%component_2_moments => null()
        this%component_1_moments => null()
    end subroutine Abstract_Mixture_Dipolar_Moments_destroy

    pure function Abstract_Mixture_Dipolar_Moments_get_num(this) result(num_moments)
        class(Abstract_Mixture_Dipolar_Moments), intent(in) :: this
        integer :: num_moments

        num_moments = this%component_1_moments%get_num() + this%component_2_moments%get_num()
    end function Abstract_Mixture_Dipolar_Moments_get_num

    pure function Abstract_Mixture_Dipolar_Moments_get(this, i_particle) result(position)
        class(Abstract_Mixture_Dipolar_Moments), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: position(num_dimensions)

        if (i_particle <= this%component_1_moments%get_num()) then
            position = this%component_1_moments%get(i_particle)
        else
            position = this%component_2_moments%get(i_particle - this%component_1_moments%get_num())
        end if
    end function Abstract_Mixture_Dipolar_Moments_get

end module class_mixture_dipolar_moments
