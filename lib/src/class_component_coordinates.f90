module class_component_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_coordinates, only: Abstract_Coordinates

implicit none

private

    type, extends(Abstract_Coordinates), abstract, public :: Abstract_Component_Coordinates
    contains
        procedure(Abstract_Component_Coordinates_set), deferred :: set
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
