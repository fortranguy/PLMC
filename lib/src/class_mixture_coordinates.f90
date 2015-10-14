module class_mixture_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use class_coordinates, only: Abstract_Coordinates
use class_component_coordinates, only: Abstract_Component_Coordinates

implicit none

private

    type, extends(Abstract_Coordinates), abstract, public :: Abstract_Mixture_Coordinates
    private
        class(Abstract_Component_Coordinates), pointer :: component_1_coordinates => null(), &
            component_2_coordinates => null()
    contains
        procedure :: construct => Abstract_Mixture_Coordinates_construct
        procedure :: destroy => Abstract_Mixture_Coordinates_destroy
        procedure :: get_num => Abstract_Mixture_Coordinates_get_num
        procedure :: get => Abstract_Mixture_Coordinates_get
    end type Abstract_Mixture_Coordinates

contains

    subroutine Abstract_Mixture_Coordinates_construct(this, component_1_coordinates, &
        component_2_coordinates)
        class(Abstract_Mixture_Coordinates), intent(out) :: this
        class(Abstract_Component_Coordinates), target, intent(in) :: component_1_coordinates, &
            component_2_coordinates

        this%component_1_coordinates => component_1_coordinates
        this%component_2_coordinates => component_2_coordinates
    end subroutine Abstract_Mixture_Coordinates_construct

    subroutine Abstract_Mixture_Coordinates_destroy(this)
        class(Abstract_Mixture_Coordinates), intent(inout) :: this

        this%component_2_coordinates => null()
        this%component_1_coordinates => null()
    end subroutine Abstract_Mixture_Coordinates_destroy

    pure function Abstract_Mixture_Coordinates_get_num(this) result(num_coordinates)
        class(Abstract_Mixture_Coordinates), intent(in) :: this
        integer :: num_coordinates

        num_coordinates = this%component_1_coordinates%get_num() + &
            this%component_2_coordinates%get_num()
    end function Abstract_Mixture_Coordinates_get_num

    pure function Abstract_Mixture_Coordinates_get(this, i_particle) result(position)
        class(Abstract_Mixture_Coordinates), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: position(num_dimensions)

        if (i_particle <= this%component_1_coordinates%get_num()) then
            position = this%component_1_coordinates%get(i_particle)
        else
            position = this%component_2_coordinates%get(i_particle - &
                this%component_1_coordinates%get_num())
        end if
    end function Abstract_Mixture_Coordinates_get

end module class_mixture_coordinates
