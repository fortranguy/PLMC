module class_random_moments

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_dipolar_spheres, only: Abstract_Dipolar_Spheres
use procedures_orientation, only: random_orientation, markov_orientation

implicit none

private

    type, abstract, public :: Abstract_Random_Moments
    private
    contains
        procedure(Abstract_Random_Moments_construct), deferred :: construct
        procedure(Abstract_Random_Moments_destroy), deferred :: destroy
        procedure(Abstract_Random_Moments_moment), deferred :: moment
        procedure(Abstract_Random_Moments_rotation), deferred :: rotation
    end type Abstract_Random_Moments

    abstract interface
    
        subroutine Abstract_Random_Moments_construct(this, dipolar_spheres, delta)
        import :: DP, Abstract_Dipolar_Spheres, Abstract_Random_Moments
            class(Abstract_Random_Moments), intent(out) :: this
            class(Abstract_Dipolar_Spheres), target, intent(in) :: dipolar_spheres
            real(DP), intent(in) :: delta
        end subroutine Abstract_Random_Moments_construct
        
        subroutine Abstract_Random_Moments_destroy(this)
        import :: Abstract_Random_Moments
            class(Abstract_Random_Moments), intent(inout) :: this
        end subroutine Abstract_Random_Moments_destroy
    
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
        procedure :: construct => Null_Random_Moments_construct
        procedure :: destroy => Null_Random_Moments_destroy
        procedure :: moment => Null_Random_Moments_moment
        procedure :: rotation => Null_Random_Moments_rotation
    end type Null_Random_Moments

    type, extends(Abstract_Random_Moments) :: Random_Moments
        real(DP) :: normed_delta
        class(Abstract_Dipolar_Spheres), pointer :: dipolar_spheres => null()
    contains
        procedure :: construct => Random_Moments_construct
        procedure :: destroy => Random_Moments_destroy
        procedure :: moment => Random_Moments_moment
        procedure :: rotation => Random_Moments_rotation
    end type Random_Moments

contains

!implementation Null_Random_Moments

    subroutine Null_Random_Moments_construct(this, dipolar_spheres, delta)
        class(Null_Random_Moments), intent(out) :: this
        class(Abstract_Dipolar_Spheres), target, intent(in) :: dipolar_spheres
        real(DP), intent(in) :: delta
        
    end subroutine Null_Random_Moments_construct
    
    subroutine Null_Random_Moments_destroy(this)
        class(Null_Random_Moments), intent(inout) :: this
        
    end subroutine Null_Random_Moments_destroy

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

    subroutine Random_Moments_construct(this, dipolar_spheres, delta)
        class(Random_Moments), intent(out) :: this
        class(Abstract_Dipolar_Spheres), target, intent(in) :: dipolar_spheres
        real(DP), intent(in) :: delta
        
        this%dipolar_spheres => dipolar_spheres
        this%normed_delta = delta / this%dipolar_spheres%get_moment_norm()
    end subroutine Random_Moments_construct
    
    subroutine Random_Moments_destroy(this)
        class(Random_Moments), intent(inout) :: this
        
        this%dipolar_spheres => null()
    end subroutine Random_Moments_destroy
    
    function Random_Moments_moment(this) result(moment)
        class(Random_Moments), intent(in) :: this
        real(DP) :: moment(num_dimensions)

        moment = this%dipolar_spheres%get_moment_norm() * random_orientation()
    end function Random_Moments_moment
    
    function Random_Moments_rotation(this, i_sphere) result(rotation)
        class(Random_Moments), intent(in) :: this
        integer, intent(in) :: i_sphere
        real(DP) :: rotation(num_dimensions)
        
        real(DP) :: orientation(num_dimensions)
        
        orientation = this%dipolar_spheres%get_moment(i_sphere)
        orientation = orientation / norm2(orientation)        
        call markov_orientation(orientation, this%normed_delta)        
        rotation = this%dipolar_spheres%get_moment_norm() * orientation
    end function Random_Moments_rotation
    
!end implementation Random_Moments

end module class_random_moments

