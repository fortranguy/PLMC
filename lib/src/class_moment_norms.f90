module class_moment_norms

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_in_range, check_positive
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
    contains
        procedure(Abstract_Moment_Norms_construct), deferred :: construct
        procedure(Abstract_Moment_Norms_destroy), deferred :: destroy
        procedure(Abstract_Moment_Norms_set), deferred :: set
        procedure(Abstract_Moment_Norms_get), deferred :: get
        procedure(Abstract_Moment_Norms_add), deferred :: add
        procedure(Abstract_Moment_Norms_remove), deferred :: remove
    end type Abstract_Moment_Norms
    
    abstract interface
    
        subroutine Abstract_Moment_Norms_construct(this, particles_num)
        import :: Abstract_Particles_Number, Abstract_Moment_Norms
            class(Abstract_Moment_Norms), intent(out) :: this
            class(Abstract_Particles_Number), target, intent(in) :: particles_num
        end subroutine Abstract_Moment_Norms_construct
        
        subroutine Abstract_Moment_Norms_destroy(this)
        import :: Abstract_Moment_Norms
            class(Abstract_Moment_Norms), intent(inout) :: this
        end subroutine Abstract_Moment_Norms_destroy
    
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
    
    type, extends(Abstract_Moment_Norms), public :: Null_Moment_Norms
    contains
        procedure :: construct => Null_Moment_Norms_construct
        procedure :: destroy => Null_Moment_Norms_destroy
        procedure :: set => Null_Moment_Norms_set
        procedure :: get => Null_Moment_Norms_get
        procedure :: add => Null_Moment_Norms_add
        procedure :: remove => Null_Moment_Norms_remove
    end type Null_Moment_Norms

    type, extends(Abstract_Moment_Norms), public :: Uniform_Moment_Norms
    private
        real(DP) :: norm
        class(Abstract_Particles_Number), pointer :: particles_num
    contains
        procedure :: construct => Uniform_Moment_Norms_construct
        procedure :: destroy => Uniform_Moment_Norms_destroy
        procedure :: set => Uniform_Moment_Norms_set
        procedure :: get => Uniform_Moment_Norms_get
        procedure :: add => Uniform_Moment_Norms_add
        procedure :: remove => Uniform_Moment_Norms_remove
    end type Uniform_Moment_Norms

contains

!implementation Null_Moment_Norms

    subroutine Null_Moment_Norms_construct(this, particles_num)
        class(Null_Moment_Norms), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_num
        
    end subroutine Null_Moment_Norms_construct
    
    subroutine Null_Moment_Norms_destroy(this)
        class(Null_Moment_Norms), intent(inout) :: this
        
    end subroutine Null_Moment_Norms_destroy

    subroutine Null_Moment_Norms_set(this, i_particle, norm)
        class(Null_Moment_Norms), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: norm
        
    end subroutine Null_Moment_Norms_set
    
    pure function Null_Moment_Norms_get(this, i_particle) result(norm)
        class(Null_Moment_Norms), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: norm
        
        norm = 0._DP
    end function Null_Moment_Norms_get
    
    subroutine Null_Moment_Norms_add(this, norm)
        class(Null_Moment_Norms), intent(inout) :: this
        real(DP), intent(in) :: norm
        
    end subroutine Null_Moment_Norms_add
    
    subroutine Null_Moment_Norms_remove(this, i_particle)
        class(Null_Moment_Norms), intent(inout) :: this
        integer, intent(in) :: i_particle
        
    end subroutine Null_Moment_Norms_remove

!end implementation Null_Moment_Norms

!implementation Uniform_Moment_Norms

    subroutine Uniform_Moment_Norms_construct(this, particles_num)
        class(Uniform_Moment_Norms), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_num
        
        this%particles_num => particles_num
    end subroutine Uniform_Moment_Norms_construct
    
    subroutine Uniform_Moment_Norms_destroy(this)
        class(Uniform_Moment_Norms), intent(inout) :: this
        
        this%particles_num => null()
    end subroutine Uniform_Moment_Norms_destroy

    subroutine Uniform_Moment_Norms_set(this, i_particle, norm)
        class(Uniform_Moment_Norms), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: norm
        
        call check_in_range("Uniform_Moment_Norms", this%particles_num%get(), &
                            "i_particle", i_particle)
        call check_positive("Uniform_Moment_Norms", "norm", norm)
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
         
        call this%set(this%particles_num%get(), norm)
    end subroutine Uniform_Moment_Norms_add
    
    subroutine Uniform_Moment_Norms_remove(this, i_particle)
         class(Uniform_Moment_Norms), intent(inout) :: this
         integer, intent(in) :: i_particle
         
        call check_in_range("Uniform_Moment_Norms", this%particles_num%get(), "i_particle", i_particle)
    end subroutine Uniform_Moment_Norms_remove
    
!end implementation Uniform_Moment_Norms

end module class_moment_norms
