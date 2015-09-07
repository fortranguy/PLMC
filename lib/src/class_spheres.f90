module class_spheres

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use module_error, only: warning_continue, error_exit
use class_particles_number, only: Abstract_Particles_Number, Abstract_Particles_Number_Pointer

implicit none

private

    type, abstract, public :: Abstract_Spheres
        type(Abstract_Particles_Number_Pointer) :: particles_num
    contains
        procedure :: construct => Abstract_Spheres_construct
        procedure :: destroy => Abstract_Spheres_destroy
        
        procedure(Abstract_Spheres_set_diameter), deferred :: set_diameter
        procedure(Abstract_Spheres_get_diameter), deferred :: get_diameter
        procedure(Abstract_Spheres_add_diameter), deferred :: add_diameter
        procedure(Abstract_Spheres_remove_diameter), deferred :: remove_diameter
    end type Abstract_Spheres
    
    type, public :: Abstract_Spheres_Pointer
        class(Abstract_Spheres), pointer :: ptr
    end type Abstract_Spheres_Pointer
    
    abstract interface
    
        subroutine Abstract_Spheres_set_diameter(this, i_particle, diameter)
        import :: DP, Abstract_Spheres
            class(Abstract_Spheres), intent(inout) :: this
            integer, intent(in) :: i_particle
            real(DP), intent(in) :: diameter
        end subroutine Abstract_Spheres_set_diameter
        
        pure function Abstract_Spheres_get_diameter(this, i_particle) result(get_diameter)
        import :: DP, Abstract_Spheres
            class(Abstract_Spheres), intent(in) :: this
            integer, intent(in) :: i_particle
            real(DP) :: get_diameter
        end function Abstract_Spheres_get_diameter
        
        subroutine Abstract_Spheres_add_diameter(this, diameter)
        import :: DP, Abstract_Spheres
            class(Abstract_Spheres), intent(inout) :: this
            real(DP), intent(in) :: diameter
        end subroutine Abstract_Spheres_add_diameter
        
        subroutine Abstract_Spheres_remove_diameter(this, i_particle)
        import :: DP, Abstract_Spheres
            class(Abstract_Spheres), intent(inout) :: this
            integer, intent(in) :: i_particle
        end subroutine Abstract_Spheres_remove_diameter
        
    end interface

    type, extends(Abstract_Spheres), public :: Uniform_Spheres
    private
        real(DP) :: diameter
    contains
        procedure :: set_diameter => Uniform_Spheres_set_diameter
        procedure :: get_diameter => Uniform_Spheres_get_diameter
        procedure :: add_diameter => Uniform_Spheres_add_diameter
        procedure :: remove_diameter => Uniform_Spheres_remove_diameter
    end type Uniform_Spheres

contains

!implementation Abstract_Spheres

    subroutine Abstract_Spheres_construct(this, particles_num)
        class(Abstract_Spheres), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_num
        
        this%particles_num%ptr => particles_num
    end subroutine Abstract_Spheres_construct
    
    subroutine Abstract_Spheres_destroy(this)
        class(Abstract_Spheres), intent(inout) :: this
        
        this%particles_num%ptr => null()
    end subroutine Abstract_Spheres_destroy

!end implementation Abstract_Spheres

!implementation Uniform_Spheres

    subroutine Uniform_Spheres_set_diameter(this, i_particle, diameter)
        class(Uniform_Spheres), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: diameter
        
        if (i_particle < 1 .or. this%particles_num%ptr%get_num() < i_particle) then
            call error_exit("Uniform_Spheres: i_particle is out of range.")
        end if
        if (diameter < 0._DP) then
            call error_exit("Uniform_Spheres: diameter is negative.")
        end if
        if (diameter < real_zero) then
            call warning_continue("Uniform_Spheres: diameter may be too small.")
        end if
        
        this%diameter = diameter
    end subroutine Uniform_Spheres_set_diameter

    pure function Uniform_Spheres_get_diameter(this, i_particle) result(get_diameter)
        class(Uniform_Spheres), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: get_diameter
        
        get_diameter = this%diameter
    end function Uniform_Spheres_get_diameter
    
    subroutine Uniform_Spheres_add_diameter(this, diameter)
         class(Uniform_Spheres), intent(inout) :: this
         real(DP), intent(in) :: diameter
         
    end subroutine Uniform_Spheres_add_diameter
    
    subroutine Uniform_Spheres_remove_diameter(this, i_particle)
         class(Uniform_Spheres), intent(inout) :: this
         integer, intent(in) :: i_particle
         
    end subroutine Uniform_Spheres_remove_diameter
    
!end implementation Uniform_Spheres

end module class_spheres
