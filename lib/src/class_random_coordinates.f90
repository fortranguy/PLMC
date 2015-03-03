module class_random_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_random_positions, only: Abstract_Random_Positions
use class_random_moments, only: Abstract_Random_Moments

implicit none

private

    type, abstract, public :: Abstract_Random_Coordinates
    private
    contains
        !Strategy Pattern
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
    
    type, extends(Abstract_Random_Coordinates), public :: Random_Coordinates
    private
        class(Abstract_Random_Positions), pointer :: rand_positions => null()
        class(Abstract_Random_Moments), pointer :: rand_moments => null()
    contains
        procedure :: position => Random_Coordinates_position
        procedure :: moment => Random_Coordinates_moment
        procedure :: move => Random_Coordinates_move
        procedure :: rotation => Random_Coordinates_rotation
    end type Random_Coordinates

contains

    function Random_Coordinates_position(this) result(position)
        class(Random_Coordinates), intent(in) :: this
        real(DP) :: position(num_dimensions)
        
        position = this%rand_positions%position()
    end function Random_Coordinates_position
    
    function Random_Coordinates_moment(this) result(moment)
        class(Random_Coordinates), intent(in) :: this
        real(DP) :: moment(num_dimensions)
        
        moment = this%rand_moments%moment()
    end function Random_Coordinates_moment
    
    function Random_Coordinates_move(this, i_sphere) result(move)
        class(Random_Coordinates), intent(in) :: this
        integer, intent(in) :: i_sphere
        real(DP) :: move(num_dimensions)
        
        move = this%rand_positions%move(i_sphere)
    end function Random_Coordinates_move
    
    function Random_Coordinates_rotation(this, i_sphere) result(rotation)
        class(Random_Coordinates), intent(in) :: this
        integer, intent(in) :: i_sphere
        real(DP) :: rotation(num_dimensions)
        
        rotation = this%rand_moments%rotation(i_sphere)
    end function Random_Coordinates_rotation

end module class_random_coordinates
