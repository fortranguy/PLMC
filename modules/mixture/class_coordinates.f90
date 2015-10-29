module class_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions

implicit none

private

    type, abstract, public :: Abstract_Coordinates
    contains
        procedure(Abstract_Coordinates_get_num), deferred :: get_num
        procedure(Abstract_Coordinates_get), deferred :: get
    end type Abstract_Coordinates

    abstract interface

        pure function Abstract_Coordinates_get_num(this) result(num_coordinates)
        import :: Abstract_Coordinates
            class(Abstract_Coordinates), intent(in) :: this
            integer :: num_coordinates
        end function Abstract_Coordinates_get_num

        pure function Abstract_Coordinates_get(this, i_particle) result(coordinate)
        import :: DP, num_dimensions, Abstract_Coordinates
            class(Abstract_Coordinates), intent(in) :: this
            integer, intent(in) :: i_particle
            real(DP) :: coordinate(num_dimensions)
        end function Abstract_Coordinates_get

    end interface

end module class_coordinates
