module class_random_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_orientation, only: gauss

implicit none

private

    type, abstract, public :: Abstract_Random_Coordinates
    private
    contains
        procedure :: position => Abstract_Random_Coordinates_position
        procedure(Abstract_Random_Coordinates_moment), deferred, nopass :: moment
        procedure :: move => Abstract_Random_Coordinates_move
        procedure(Abstract_Random_Coordinates_rotation), deferred :: rotation
    end type Abstract_Random_Coordinates

    abstract interface
    
        function Abstract_Random_Coordinates_moment(this) result(moment)
        import :: DP, num_dimensions, Abstract_Random_Coordinates
            class(Abstract_Random_Coordinates), intent(in) :: this
            real(DP) :: moment(num_dimensions)
        end function Abstract_Random_Coordinates_moment
        
        function Abstract_Random_Coordinates_rotation(this, i_sphere) result(rotation)
        import :: DP, num_dimensions
            class(Abstract_Random_Coordinates), intent(in) :: this
            real(DP) :: rotation
        end function Abstract_Random_Coordinates_rotation
        
    end interface

    type, extends(Abstract_Random_Coordinates) :: Random_Positions
    contains
        procedure :: moment => Random_Positions_moment
        procedure :: rotation => Random_Positions_rotation
    end type Random_Positions

    type, extends(Abstract_Random_Coordinates) :: Random_Coordinates
    contains
        procedure :: moment => Random_Coordinates_moment
        procedure :: rotation => Random_Coordinates_rotation
    end type Random_Coordinates

contains

!implementation Abstract_Random_Coordinates

    pure function Abstract_Random_Coordinates_position(this)
        class(Abstract_Random_Coordinates), intent(in) :: this
    end function Abstract_Random_Coordinates_position

!end implementation Abstract_Random_Coordinates

!implementation Random_Positions
    function Random_Positions_moment() result(moment)
        real(DP) :: moment(num_dimensions)
        
        moment = 0._DP
    end function Random_Positions_moment
!end implementation Random_Positions

!implementation Random_Coordinates
    !> From SMAC, Algorithm 1.23, p. 43
    function Random_Coordinates_moment() result(moment)
        real(DP) :: moment(num_dimensions)

        integer :: i_dimension

        do i_dimension = 1, num_dimensions
            moment(i_dimension) = gauss()
        end do
        moment = moment / norm2(moment)
    end function Random_Coordinates_moment
!end implementation Random_Coordinates

end module class_random_coordinates
