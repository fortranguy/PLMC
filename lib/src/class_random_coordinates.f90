module class_random_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_orientation, only: gauss

implicit none

private

    type, abstract, public :: Abstract_Random_Coordinates
    private
    contains
        procedure(Abstract_Random_Coordinates_position), deferred :: position
        procedure(Abstract_Random_Coordinates_moment), deferred :: moment
        procedure(Abstract_Random_Coordinates_move), deferred :: move
        procedure(Abstract_Random_Coordinates_rotation), deferred :: rotation
    end type Abstract_Random_Coordinates

    abstract interface
    
        function Abstract_Random_Coordinates_position(this) result(position)
        import :: DP, num_dimensions, Abstract_Random_Coordinates
            class(Abstract_Random_Coordinates), intent(in) :: this
            real(DP) :: position(num_dimensions)
        end function Abstract_Random_Coordinates_position
    
        function Abstract_Random_Coordinates_moment(this) result(moment)
        import :: DP, num_dimensions, Abstract_Random_Coordinates
            class(Abstract_Random_Coordinates), intent(in) :: this
            real(DP) :: moment(num_dimensions)
        end function Abstract_Random_Coordinates_moment
        
        function Abstract_Random_Coordinates_move(this, i_sphere) result(move)
        import :: DP, num_dimensions, Abstract_Random_Coordinates
            class(Abstract_Random_Coordinates), intent(in) :: this
            integer, intent(in) :: i_sphere
            real(DP) :: move(num_dimensions)
        end function Abstract_Random_Coordinates_move
        
        function Abstract_Random_Coordinates_rotation(this, i_sphere) result(rotation)
        import :: DP, num_dimensions, Abstract_Random_Coordinates
            class(Abstract_Random_Coordinates), intent(in) :: this
            integer, intent(in) :: i_sphere
            real(DP) :: rotation(num_dimensions)
        end function Abstract_Random_Coordinates_rotation
        
    end interface

contains

end module class_random_coordinates
