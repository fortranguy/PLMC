module class_random_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_random_positions, only: Abstract_Random_Positions
use class_random_moments, only: Abstract_Random_Moments

implicit none

private

    !> Strategy Pattern
    type, abstract, public :: Abstract_Random_Coordinates
    private
        class(Abstract_Random_Positions), pointer :: rand_positions => null()
        class(Abstract_Random_Moments), pointer :: rand_moments => null()
    contains
        procedure :: construct => Abstract_Random_Coordinates_construct
        procedure :: destroy => Abstract_Random_Coordinates_destroy

        procedure :: position => Abstract_Random_Coordinates_position
        procedure :: moment => Abstract_Random_Coordinates_moment
        procedure :: move => Abstract_Random_Coordinates_move
        procedure :: rotation => Abstract_Random_Coordinates_rotation
    end type Abstract_Random_Coordinates
    
    type, extends(Abstract_Random_Coordinates), public :: Random_Coordinates
    
    end type Random_Coordinates

contains

    subroutine Abstract_Random_Coordinates_construct(this, rand_positions, rand_moments)
        class(Abstract_Random_Coordinates), intent(out) :: this
        class(Abstract_Random_Positions), target, intent(in) :: rand_positions
        class(Abstract_Random_Moments), target, intent(in) :: rand_moments
        
        this%rand_positions => rand_positions
        this%rand_moments => rand_moments
    end subroutine Abstract_Random_Coordinates_construct
    
    subroutine Abstract_Random_Coordinates_destroy(this)
        class(Abstract_Random_Coordinates), intent(inout) :: this
        
        this%rand_moments => null()
        this%rand_positions => null()
    end subroutine Abstract_Random_Coordinates_destroy

    function Abstract_Random_Coordinates_position(this) result(position)
        class(Abstract_Random_Coordinates), intent(in) :: this
        real(DP) :: position(num_dimensions)
        
        position = this%rand_positions%position()
    end function Abstract_Random_Coordinates_position
    
    function Abstract_Random_Coordinates_moment(this) result(moment)
        class(Abstract_Random_Coordinates), intent(in) :: this
        real(DP) :: moment(num_dimensions)
        
        moment = this%rand_moments%moment()
    end function Abstract_Random_Coordinates_moment
    
    function Abstract_Random_Coordinates_move(this, i_sphere) result(move)
        class(Abstract_Random_Coordinates), intent(in) :: this
        integer, intent(in) :: i_sphere
        real(DP) :: move(num_dimensions)
        
        move = this%rand_positions%move(i_sphere)
    end function Abstract_Random_Coordinates_move
    
    function Abstract_Random_Coordinates_rotation(this, i_sphere) result(rotation)
        class(Abstract_Random_Coordinates), intent(in) :: this
        integer, intent(in) :: i_sphere
        real(DP) :: rotation(num_dimensions)
        
        rotation = this%rand_moments%rotation(i_sphere)
    end function Abstract_Random_Coordinates_rotation

end module class_random_coordinates
