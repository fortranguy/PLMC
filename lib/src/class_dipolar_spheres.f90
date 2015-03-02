module class_dipolar_spheres

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions

implicit none

private

    type, abstract, public :: Abstract_Dipolar_Spheres
    private
        character(len=:), allocatable :: name
        real(DP) :: diameter
        integer ::  num
        real(DP), allocatable :: positions(:, :)
    contains
        procedure, non_overridable :: get_name => Abstract_Dipolar_Spheres_get_name
        procedure, non_overridable :: get_diameter => Abstract_Dipolar_Spheres_get_diameter
        procedure, non_overridable :: get_num => Abstract_Dipolar_Spheres_get_num
        procedure, non_overridable :: get_position => Abstract_Dipolar_Spheres_get_position
        procedure, non_overridable :: set_position => Abstract_Dipolar_Spheres_set_position
        procedure(Abstract_Dipolar_Spheres_get_moment_norm), deferred :: get_moment_norm
        procedure(Abstract_Dipolar_Spheres_get_moment), deferred :: get_moment
        procedure(Abstract_Dipolar_Spheres_set_moment), deferred :: set_moment
    end type Abstract_Dipolar_Spheres
    
    abstract interface
    
        pure function Abstract_Dipolar_Spheres_get_moment_norm(this) result(get_moment_norm)
        import :: Abstract_Dipolar_Spheres, DP
            class(Abstract_Dipolar_Spheres), intent(in) :: this
            real(DP) :: get_moment_norm
        end function Abstract_Dipolar_Spheres_get_moment_norm
        
        pure function Abstract_Dipolar_Spheres_get_moment(this) result(get_moment)
        import :: DP, num_dimensions, Abstract_Dipolar_Spheres
            class(Abstract_Dipolar_Spheres), intent(in) :: this
            real(DP) :: get_moment(num_dimensions)
        end function Abstract_Dipolar_Spheres_get_moment
        
        pure subroutine Abstract_Dipolar_Spheres_set_moment(this, i_sphere, moment)
            import :: DP, num_dimensions, Abstract_Dipolar_Spheres
            class(Abstract_Dipolar_Spheres), intent(inout) :: this
            integer, intent(in) :: i_sphere
            real(DP), intent(in) :: moment(num_dimensions)
        end subroutine Abstract_Dipolar_Spheres_set_moment
    
    end interface
    
contains

    pure function Abstract_Dipolar_Spheres_get_name(this) result(get_name)
        class(Abstract_Dipolar_Spheres), intent(in) :: this
        character(len=:), allocatable :: get_name
        
        get_name = this%name
    end function Abstract_Dipolar_Spheres_get_name
    
    pure function Abstract_Dipolar_Spheres_get_diameter(this) result(get_diameter)
        class(Abstract_Dipolar_Spheres), intent(in) :: this
        real(DP) :: get_diameter
        
        get_diameter = this%diameter
    end function Abstract_Dipolar_Spheres_get_diameter
    
    pure function Abstract_Dipolar_Spheres_get_num(this) result(get_num)
        class(Abstract_Dipolar_Spheres), intent(in) :: this
        integer :: get_num
        
        get_num = this%num
    end function Abstract_Dipolar_Spheres_get_num
    
    pure function Abstract_Dipolar_Spheres_get_position(this, i_sphere) result(get_position)
        class(Abstract_Dipolar_Spheres), intent(in) :: this
        integer, intent(in) :: i_sphere
        real(DP) :: get_position(num_dimensions)
        
        get_position = this%positions(:, i_sphere)
    end function Abstract_Dipolar_Spheres_get_position
    
    pure subroutine Abstract_Dipolar_Spheres_set_position(this, i_sphere, position)
        class(Abstract_Dipolar_Spheres), intent(inout) :: this
        integer, intent(in) :: i_sphere
        real(DP), intent(in) :: position(num_dimensions)
        
        this%positions(:, i_sphere) = position
    end subroutine Abstract_Dipolar_Spheres_set_position

end module class_dipolar_spheres
