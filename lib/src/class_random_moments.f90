module class_random_moments

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_orientation, only: gauss
use class_dipolar_spheres, only: Abstract_Dipolar_Spheres

implicit none

private

    type, abstract, public :: Abstract_Random_Moments
    private
    contains
        procedure(Abstract_Random_Moments_moment), deferred :: moment
        procedure(Abstract_Random_Moments_rotation), deferred :: rotation
    end type Abstract_Random_Moments

    abstract interface
    
        function Abstract_Random_Moments_moment(this) result(moment)
        import :: DP, num_dimensions, Abstract_Random_Moments
            class(Abstract_Random_Moments), intent(in) :: this
            real(DP) :: moment(num_dimensions)
        end function Abstract_Random_Moments_moment
        
        function Abstract_Random_Moments_rotation(this, i_sphere) result(rotation)
        import :: DP, num_dimensions, Abstract_Random_Moments
            class(Abstract_Random_Moments), intent(in) :: this
            integer, intent(in) :: i_sphere
            real(DP) :: rotation(num_dimensions)
        end function Abstract_Random_Moments_rotation
           
    end interface

    type, extends(Abstract_Random_Moments) :: Null_Random_Moments
    contains
        procedure :: moment => Null_Random_Moments_moment
        procedure :: rotation => Null_Random_Moments_rotation
    end type Null_Random_Moments

    type, extends(Abstract_Random_Moments) :: Random_Moments
        class(Abstract_Dipolar_Spheres), pointer :: dipolar_spheres => null()
    contains
        procedure :: moment => Random_Moments_moment
        procedure :: rotation => Random_Moments_rotation
    end type Random_Moments

contains

!implementation Null_Random_Moments

    function Null_Random_Moments_moment(this) result(moment)
        class(Null_Random_Moments), intent(in) :: this
        real(DP) :: moment(num_dimensions)
        
        moment = 0._DP
    end function Null_Random_Moments_moment
    
    function Null_Random_Moments_rotation(this, i_sphere) result(rotation)
        class(Null_Random_Moments), intent(in) :: this
        integer, intent(in) :: i_sphere
        real(DP) :: rotation(num_dimensions)
        
        rotation = 0._DP
    end function Null_Random_Moments_rotation
    
!end implementation Null_Random_Moments

!implementation Random_Moments

    !> From SMAC, Algorithm 1.23, p. 43
    function Random_Moments_moment(this) result(moment)
        class(Random_Moments), intent(in) :: this
        real(DP) :: moment(num_dimensions)

        integer :: i_dimension

        do i_dimension = 1, num_dimensions
            moment(i_dimension) = gauss()
        end do
        moment = moment / norm2(moment) ! * moment_norm()
    end function Random_Moments_moment
    
    function Random_Moments_rotation(this, i_sphere) result(rotation)
        class(Random_Moments), intent(in) :: this
        integer, intent(in) :: i_sphere
        real(DP) :: rotation(num_dimensions)
    end function Random_Moments_rotation
    
!end implementation Random_Moments

end module class_random_moments
