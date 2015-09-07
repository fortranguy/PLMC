module class_dipolar_moment_norms

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use module_error, only: warning_continue, error_exit
use class_particles_number, only: Abstract_Particles_Number, Abstract_Particles_Number_ptr

implicit none

private

    type, abstract, public :: Abstract_Dipolar_Moment_Norms
        type(Abstract_Particles_Number_ptr) :: particles_num
    contains
        procedure :: construct => Abstract_Dipolar_Moment_Norms_construct
        procedure :: destroy => Abstract_Dipolar_Moment_Norms_destroy
        
        procedure(Abstract_Dipolar_Moment_Norms_set_norm), deferred :: set_norm
        procedure(Abstract_Dipolar_Moment_Norms_get_norm), deferred :: get_norm
        procedure(Abstract_Dipolar_Moment_Norms_add_norm), deferred :: add_norm
        procedure(Abstract_Dipolar_Moment_Norms_remove_norm), deferred :: remove_norm
    end type Abstract_Dipolar_Moment_Norms
    
    abstract interface
    
        subroutine Abstract_Dipolar_Moment_Norms_set_norm(this, i_particle, norm)
        import :: DP, Abstract_Dipolar_Moment_Norms
            class(Abstract_Dipolar_Moment_Norms), intent(inout) :: this
            integer, intent(in) :: i_particle
            real(DP), intent(in) :: norm
        end subroutine Abstract_Dipolar_Moment_Norms_set_norm
        
        pure function Abstract_Dipolar_Moment_Norms_get_norm(this, i_particle) result(get_norm)
        import :: DP, Abstract_Dipolar_Moment_Norms
            class(Abstract_Dipolar_Moment_Norms), intent(in) :: this
            integer, intent(in) :: i_particle
            real(DP) :: get_norm
        end function Abstract_Dipolar_Moment_Norms_get_norm
        
        subroutine Abstract_Dipolar_Moment_Norms_add_norm(this, norm)
        import :: DP, Abstract_Dipolar_Moment_Norms
            class(Abstract_Dipolar_Moment_Norms), intent(inout) :: this
            real(DP), intent(in) :: norm
        end subroutine Abstract_Dipolar_Moment_Norms_add_norm
        
        subroutine Abstract_Dipolar_Moment_Norms_remove_norm(this, i_particle)
        import :: DP, Abstract_Dipolar_Moment_Norms
            class(Abstract_Dipolar_Moment_Norms), intent(inout) :: this
            integer, intent(in) :: i_particle
        end subroutine Abstract_Dipolar_Moment_Norms_remove_norm
        
    end interface

    type, extends(Abstract_Dipolar_Moment_Norms), public :: Uniform_Dipolar_Moment_Norms
    private
        real(DP) :: norm
    contains
        procedure :: set_norm => Uniform_Dipolar_Moment_Norms_set_norm
        procedure :: get_norm => Uniform_Dipolar_Moment_Norms_get_norm
        procedure :: add_norm => Uniform_Dipolar_Moment_Norms_add_norm
        procedure :: remove_norm => Uniform_Dipolar_Moment_Norms_remove_norm
    end type Uniform_Dipolar_Moment_Norms

contains

!implementation Abstract_Dipolar_Moment_Norms

    subroutine Abstract_Dipolar_Moment_Norms_construct(this, particles_num)
        class(Abstract_Dipolar_Moment_Norms), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_num
        
        this%particles_num%ptr => particles_num
    end subroutine Abstract_Dipolar_Moment_Norms_construct
    
    subroutine Abstract_Dipolar_Moment_Norms_destroy(this)
        class(Abstract_Dipolar_Moment_Norms), intent(inout) :: this
        
        this%particles_num%ptr => null()
    end subroutine Abstract_Dipolar_Moment_Norms_destroy

!end implementation Abstract_Dipolar_Moment_Norms

!implementation Uniform_Dipolar_Moment_Norms

    subroutine Uniform_Dipolar_Moment_Norms_set_norm(this, i_particle, norm)
        class(Uniform_Dipolar_Moment_Norms), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: norm
        
        if (i_particle < 1 .or. this%particles_num%ptr%get_num() < i_particle) then
            call error_exit("Uniform_Dipolar_Moment_Norms: i_particle is out of range.")
        end if
        if (norm < 0._DP) then
            call error_exit("Uniform_Dipolar_Moment_Norms: norm is negative.")
        end if
        if (norm < real_zero) then
            call warning_continue("Uniform_Dipolar_Moment_Norms: norm may be too small.")
        end if
        
        this%norm = norm
    end subroutine Uniform_Dipolar_Moment_Norms_set_norm

    pure function Uniform_Dipolar_Moment_Norms_get_norm(this, i_particle) result(get_norm)
        class(Uniform_Dipolar_Moment_Norms), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: get_norm
        
        get_norm = this%norm
    end function Uniform_Dipolar_Moment_Norms_get_norm
    
    subroutine Uniform_Dipolar_Moment_Norms_add_norm(this, norm)
         class(Uniform_Dipolar_Moment_Norms), intent(inout) :: this
         real(DP), intent(in) :: norm
         
    end subroutine Uniform_Dipolar_Moment_Norms_add_norm
    
    subroutine Uniform_Dipolar_Moment_Norms_remove_norm(this, i_particle)
         class(Uniform_Dipolar_Moment_Norms), intent(inout) :: this
         integer, intent(in) :: i_particle
         
    end subroutine Uniform_Dipolar_Moment_Norms_remove_norm
    
!end implementation Uniform_Dipolar_Moment_Norms

end module class_dipolar_moment_norms
