module class_mixture_positions

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use class_coordinates, only: Abstract_Coordinates
use class_component_positions, only: Abstract_Component_Positions

implicit none

private

    type, extends(Abstract_Coordinates), abstract, public :: Abstract_Mixture_Positions
    private
        class(Abstract_Component_Positions), pointer :: component_1_positions => null(), &
            component_2_positions => null()
    contains
        procedure :: construct => Abstract_Mixture_Positions_construct
        procedure :: destroy => Abstract_Mixture_Positions_destroy
        procedure :: get_num => Abstract_Mixture_Positions_get_num
        procedure :: get => Abstract_Mixture_Positions_get
    end type Abstract_Mixture_Positions

contains

    subroutine Abstract_Mixture_Positions_construct(this, component_1_positions, &
        component_2_positions)
        class(Abstract_Mixture_Positions), intent(out) :: this
        class(Abstract_Component_Positions), target, intent(in) :: component_1_positions, &
            component_2_positions

        this%component_1_positions => component_1_positions
        this%component_2_positions => component_2_positions
    end subroutine Abstract_Mixture_Positions_construct

    subroutine Abstract_Mixture_Positions_destroy(this)
        class(Abstract_Mixture_Positions), intent(inout) :: this

        this%component_2_positions => null()
        this%component_1_positions => null()
    end subroutine Abstract_Mixture_Positions_destroy

    pure function Abstract_Mixture_Positions_get_num(this) result(num_positions)
        class(Abstract_Mixture_Positions), intent(in) :: this
        integer :: num_positions

        num_positions = this%component_1_positions%get_num() + this%component_2_positions%get_num()
    end function Abstract_Mixture_Positions_get_num

    pure function Abstract_Mixture_Positions_get(this, i_particle) result(position)
        class(Abstract_Mixture_Positions), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: position(num_dimensions)

        if (i_particle <= this%component_1_positions%get_num()) then
            position = this%component_1_positions%get(i_particle)
        else
            position = this%component_2_positions%get(i_particle - &
                this%component_1_positions%get_num())
        end if
    end function Abstract_Mixture_Positions_get

end module class_mixture_positions
