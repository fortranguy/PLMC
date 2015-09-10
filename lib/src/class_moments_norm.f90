module class_moments_norm

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use procedures_errors, only: warning_continue
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

    type, abstract, public :: Abstract_Moments_Norm
    contains
        procedure(Abstract_Moments_Norm_construct), deferred :: construct
        procedure(Abstract_Moments_Norm_destroy), deferred :: destroy
        procedure(Abstract_Moments_Norm_set), deferred :: set
        procedure(Abstract_Moments_Norm_get), deferred :: get
        procedure(Abstract_Moments_Norm_add), deferred :: add
        procedure(Abstract_Moments_Norm_remove), deferred :: remove
    end type Abstract_Moments_Norm
    
    abstract interface
    
        subroutine Abstract_Moments_Norm_construct(this, particles_number)
        import :: Abstract_Particles_Number, Abstract_Moments_Norm
            class(Abstract_Moments_Norm), intent(out) :: this
            class(Abstract_Particles_Number), target, intent(in) :: particles_number
        end subroutine Abstract_Moments_Norm_construct
        
        subroutine Abstract_Moments_Norm_destroy(this)
        import :: Abstract_Moments_Norm
            class(Abstract_Moments_Norm), intent(inout) :: this
        end subroutine Abstract_Moments_Norm_destroy
    
        subroutine Abstract_Moments_Norm_set(this, i_particle, norm)
        import :: DP, Abstract_Moments_Norm
            class(Abstract_Moments_Norm), intent(inout) :: this
            integer, intent(in) :: i_particle
            real(DP), intent(in) :: norm
        end subroutine Abstract_Moments_Norm_set
        
        pure function Abstract_Moments_Norm_get(this, i_particle) result(norm)
        import :: DP, Abstract_Moments_Norm
            class(Abstract_Moments_Norm), intent(in) :: this
            integer, intent(in) :: i_particle
            real(DP) :: norm
        end function Abstract_Moments_Norm_get
        
        subroutine Abstract_Moments_Norm_add(this, norm)
        import :: DP, Abstract_Moments_Norm
            class(Abstract_Moments_Norm), intent(inout) :: this
            real(DP), intent(in) :: norm
        end subroutine Abstract_Moments_Norm_add
        
        subroutine Abstract_Moments_Norm_remove(this, i_particle)
        import :: DP, Abstract_Moments_Norm
            class(Abstract_Moments_Norm), intent(inout) :: this
            integer, intent(in) :: i_particle
        end subroutine Abstract_Moments_Norm_remove
        
    end interface
    
    type, extends(Abstract_Moments_Norm), public :: Null_Moments_Norm
    contains
        procedure :: construct => Null_Moments_Norm_construct
        procedure :: destroy => Null_Moments_Norm_destroy
        procedure :: set => Null_Moments_Norm_set
        procedure :: get => Null_Moments_Norm_get
        procedure :: add => Null_Moments_Norm_add
        procedure :: remove => Null_Moments_Norm_remove
    end type Null_Moments_Norm

    type, extends(Abstract_Moments_Norm), public :: Uniform_Moments_Norm
    private
        real(DP) :: norm
        logical :: is_set = .false.
        class(Abstract_Particles_Number), pointer :: particles_number
    contains
        procedure :: construct => Uniform_Moments_Norm_construct
        procedure :: destroy => Uniform_Moments_Norm_destroy
        procedure :: set => Uniform_Moments_Norm_set
        procedure :: get => Uniform_Moments_Norm_get
        procedure :: add => Uniform_Moments_Norm_add
        procedure :: remove => Uniform_Moments_Norm_remove
    end type Uniform_Moments_Norm

contains

!implementation Null_Moments_Norm

    subroutine Null_Moments_Norm_construct(this, particles_number)
        class(Null_Moments_Norm), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_number
        
    end subroutine Null_Moments_Norm_construct
    
    subroutine Null_Moments_Norm_destroy(this)
        class(Null_Moments_Norm), intent(inout) :: this
        
    end subroutine Null_Moments_Norm_destroy

    subroutine Null_Moments_Norm_set(this, i_particle, norm)
        class(Null_Moments_Norm), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: norm
        
    end subroutine Null_Moments_Norm_set
    
    pure function Null_Moments_Norm_get(this, i_particle) result(norm)
        class(Null_Moments_Norm), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: norm
        
        norm = 0._DP
    end function Null_Moments_Norm_get
    
    subroutine Null_Moments_Norm_add(this, norm)
        class(Null_Moments_Norm), intent(inout) :: this
        real(DP), intent(in) :: norm
        
    end subroutine Null_Moments_Norm_add
    
    subroutine Null_Moments_Norm_remove(this, i_particle)
        class(Null_Moments_Norm), intent(inout) :: this
        integer, intent(in) :: i_particle
        
    end subroutine Null_Moments_Norm_remove

!end implementation Null_Moments_Norm

!implementation Uniform_Moments_Norm

    subroutine Uniform_Moments_Norm_construct(this, particles_number)
        class(Uniform_Moments_Norm), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_number
        
        this%particles_number => particles_number
    end subroutine Uniform_Moments_Norm_construct
    
    subroutine Uniform_Moments_Norm_destroy(this)
        class(Uniform_Moments_Norm), intent(inout) :: this
        
        this%particles_number => null()
    end subroutine Uniform_Moments_Norm_destroy

    subroutine Uniform_Moments_Norm_set(this, i_particle, norm)
        class(Uniform_Moments_Norm), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: norm
        
        call check_in_range("Uniform_Moments_Norm", this%particles_number%get(), &
                            "i_particle", i_particle)
        call check_positive("Uniform_Moments_Norm", "norm", norm)
        if (.not. this%is_set) then
            this%is_set = .true.
        else if (abs(norm - this%norm) > real_zero) then
            call warning_continue("Uniform_Moments_Norm: setting norm is different.")
        end if
        this%norm = norm
    end subroutine Uniform_Moments_Norm_set

    pure function Uniform_Moments_Norm_get(this, i_particle) result(norm)
        class(Uniform_Moments_Norm), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: norm
        
        norm = this%norm
    end function Uniform_Moments_Norm_get

    subroutine Uniform_Moments_Norm_add(this, norm)
        class(Uniform_Moments_Norm), intent(inout) :: this
        real(DP), intent(in) :: norm
        
        if (abs(norm - this%norm) > real_zero) then
            call warning_continue("Uniform_Moments_Norm: adding norm is different.")
        end if
        call this%set(this%particles_number%get(), norm)
    end subroutine Uniform_Moments_Norm_add
    
    subroutine Uniform_Moments_Norm_remove(this, i_particle)
         class(Uniform_Moments_Norm), intent(inout) :: this
         integer, intent(in) :: i_particle
         
        call check_in_range("Uniform_Moments_Norm", this%particles_number%get(), &
                            "i_particle", i_particle)
    end subroutine Uniform_Moments_Norm_remove
    
!end implementation Uniform_Moments_Norm

end module class_moments_norm
