module class_moment_norms

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive_real
use class_particles_number, only: Abstract_Particles_Number

implicit none

private

    !> A hack to count the number of objects
    !num_dip_spheres_objects = num_dip_spheres_objects + 1
    !if (num_dip_spheres_objects == 1) then
    !    if (abs(this%moment_norm - 1._DP) > real_zero) then
    !        write(error_unit, *) data_field, " must be 1.0 since it is the unit of moment."
    !        error stop
    !    end if
    !end if

    type, abstract, public :: Abstract_Moment_Norms
        class(Abstract_Particles_Number), pointer :: particles_num
    contains
        procedure :: construct => Abstract_Moment_Norms_construct
        procedure :: destroy => Abstract_Moment_Norms_destroy
        
        procedure(Abstract_Moment_Norms_set), deferred :: set
        procedure(Abstract_Moment_Norms_get), deferred :: get
        procedure(Abstract_Moment_Norms_add), deferred :: add
        procedure(Abstract_Moment_Norms_remove), deferred :: remove
    end type Abstract_Moment_Norms
    
    abstract interface
    
        subroutine Abstract_Moment_Norms_set(this, i_particle, norm)
        import :: DP, Abstract_Moment_Norms
            class(Abstract_Moment_Norms), intent(inout) :: this
            integer, intent(in) :: i_particle
            real(DP), intent(in) :: norm
        end subroutine Abstract_Moment_Norms_set
        
        pure function Abstract_Moment_Norms_get(this, i_particle) result(norm)
        import :: DP, Abstract_Moment_Norms
            class(Abstract_Moment_Norms), intent(in) :: this
            integer, intent(in) :: i_particle
            real(DP) :: norm
        end function Abstract_Moment_Norms_get
        
        subroutine Abstract_Moment_Norms_add(this, norm)
        import :: DP, Abstract_Moment_Norms
            class(Abstract_Moment_Norms), intent(inout) :: this
            real(DP), intent(in) :: norm
        end subroutine Abstract_Moment_Norms_add
        
        subroutine Abstract_Moment_Norms_remove(this, i_particle)
        import :: DP, Abstract_Moment_Norms
            class(Abstract_Moment_Norms), intent(inout) :: this
            integer, intent(in) :: i_particle
        end subroutine Abstract_Moment_Norms_remove
        
    end interface

    type, extends(Abstract_Moment_Norms), public :: Uniform_Moment_Norms
    private
        real(DP) :: norm
    contains
        procedure :: set => Uniform_Moment_Norms_set
        procedure :: get => Uniform_Moment_Norms_get
        procedure :: add => Uniform_Moment_Norms_add
        procedure :: remove => Uniform_Moment_Norms_remove
    end type Uniform_Moment_Norms

contains

!implementation Abstract_Moment_Norms

    subroutine Abstract_Moment_Norms_construct(this, particles_num)
        class(Abstract_Moment_Norms), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_num
        
        this%particles_num => particles_num
    end subroutine Abstract_Moment_Norms_construct
    
    subroutine Abstract_Moment_Norms_destroy(this)
        class(Abstract_Moment_Norms), intent(inout) :: this
        
        this%particles_num => null()
    end subroutine Abstract_Moment_Norms_destroy

!end implementation Abstract_Moment_Norms

!implementation Uniform_Moment_Norms

    subroutine Uniform_Moment_Norms_set(this, i_particle, norm)
        class(Uniform_Moment_Norms), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: norm
        
        if (i_particle < 1 .or. this%particles_num%get() < i_particle) then
            call error_exit("Uniform_Moment_Norms: i_particle is out of range.")
        end if
        call check_positive_real("Uniform_Moment_Norms", "norm", norm)
        this%norm = norm
    end subroutine Uniform_Moment_Norms_set

    pure function Uniform_Moment_Norms_get(this, i_particle) result(norm)
        class(Uniform_Moment_Norms), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: norm
        
        norm = this%norm
    end function Uniform_Moment_Norms_get
    
    subroutine Uniform_Moment_Norms_add(this, norm)
         class(Uniform_Moment_Norms), intent(inout) :: this
         real(DP), intent(in) :: norm
         
    end subroutine Uniform_Moment_Norms_add
    
    subroutine Uniform_Moment_Norms_remove(this, i_particle)
         class(Uniform_Moment_Norms), intent(inout) :: this
         integer, intent(in) :: i_particle
         
    end subroutine Uniform_Moment_Norms_remove
    
!end implementation Uniform_Moment_Norms

end module class_moment_norms
