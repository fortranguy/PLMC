module class_component_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions

implicit none

private

    type, abstract, public :: Abstract_Component_Coordinates
    contains
        procedure(Abstract_Component_Coordinates_set), deferred :: set
        procedure(Abstract_Component_Coordinates_get_num), deferred :: get_num
        procedure(Abstract_Component_Coordinates_get),deferred :: get
        procedure(Abstract_Component_Coordinates_add), deferred :: add
        procedure(Abstract_Component_Coordinates_remove), deferred :: remove
    end type Abstract_Component_Coordinates

    abstract interface

        subroutine Abstract_Component_Coordinates_set(this, i_particle, vector)
        import :: DP, Abstract_Component_Coordinates
            class(Abstract_Component_Coordinates), intent(inout) :: this
            integer, intent(in) :: i_particle
            real(DP), intent(in) :: vector(:)
        end subroutine Abstract_Component_Coordinates_set

        pure function Abstract_Component_Coordinates_get_num(this) result(num_vectors)
        import :: Abstract_Component_Coordinates
            class(Abstract_Component_Coordinates), intent(in) :: this
            integer :: num_vectors
        end function Abstract_Component_Coordinates_get_num

        pure function Abstract_Component_Coordinates_get(this, i_particle) result(vector)
        import :: DP, num_dimensions, Abstract_Component_Coordinates
            class(Abstract_Component_Coordinates), intent(in) :: this
            integer, intent(in) :: i_particle
            real(DP) :: vector(num_dimensions)
        end function Abstract_Component_Coordinates_get

        subroutine Abstract_Component_Coordinates_add(this, vector)
        import :: DP, Abstract_Component_Coordinates
            class(Abstract_Component_Coordinates), intent(inout) :: this
            real(DP), intent(in) :: vector(:)
        end subroutine Abstract_Component_Coordinates_add

        subroutine Abstract_Component_Coordinates_remove(this, i_particle)
        import :: Abstract_Component_Coordinates
            class(Abstract_Component_Coordinates), intent(inout) :: this
            integer, intent(in) :: i_particle
        end subroutine Abstract_Component_Coordinates_remove

    end interface

end module class_component_coordinates
