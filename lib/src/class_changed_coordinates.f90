module class_changed_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
implicit none

private

    type, abstract, public :: Abstract_Changed_Coordinates
    contains
        procedure(Abstract_Changed_Coordinates_increase_delta), deferred :: increase_delta
        procedure(Abstract_Changed_Coordinates_decrease_delta), deferred :: decrease_delta
        procedure(Abstract_Changed_Coordinates_get), deferred :: get
    end type

    abstract interface

        subroutine Abstract_Changed_Coordinates_increase_delta(this)
        import :: Abstract_Changed_Coordinates
            class(Abstract_Changed_Coordinates), intent(inout) :: this
        end subroutine Abstract_Changed_Coordinates_increase_delta

        subroutine Abstract_Changed_Coordinates_decrease_delta(this)
        import :: Abstract_Changed_Coordinates
            class(Abstract_Changed_Coordinates), intent(inout) :: this
        end subroutine Abstract_Changed_Coordinates_decrease_delta

        function Abstract_Changed_Coordinates_get(this, i_particle) result(changed)
            import :: DP, num_dimensions, Abstract_Changed_Coordinates
            real(DP) :: changed(num_dimensions)
            class(Abstract_Changed_Coordinates), intent(in) :: this
            integer, intent(in) :: i_particle
        end function Abstract_Changed_Coordinates_get

    end interface

end module class_changed_coordinates
