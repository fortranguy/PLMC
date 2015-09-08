module class_diameters

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use module_error, only: warning_continue, error_exit
use class_particles_number, only: Abstract_Particles_Number, Abstract_Particles_Number_Pointer

implicit none

private

    type, abstract, public :: Abstract_Diameters
        type(Abstract_Particles_Number_Pointer) :: particles_num
    contains
        procedure :: construct => Abstract_Diameters_construct
        procedure :: destroy => Abstract_Diameters_destroy
        
        procedure(Abstract_Diameters_set), deferred :: set
        procedure(Abstract_Diameters_get), deferred :: get
        procedure(Abstract_Diameters_add), deferred :: add
        procedure(Abstract_Diameters_remove), deferred :: remove
    end type Abstract_Diameters
    
    type, public :: Abstract_Diameters_Pointer
        class(Abstract_Diameters), pointer :: ptr
    end type Abstract_Diameters_Pointer
    
    abstract interface
    
        subroutine Abstract_Diameters_set(this, i_particle, diameter)
        import :: DP, Abstract_Diameters
            class(Abstract_Diameters), intent(inout) :: this
            integer, intent(in) :: i_particle
            real(DP), intent(in) :: diameter
        end subroutine Abstract_Diameters_set
        
        pure function Abstract_Diameters_get(this, i_particle) result(diameter)
        import :: DP, Abstract_Diameters
            class(Abstract_Diameters), intent(in) :: this
            integer, intent(in) :: i_particle
            real(DP) :: diameter
        end function Abstract_Diameters_get
        
        subroutine Abstract_Diameters_add(this, diameter)
        import :: DP, Abstract_Diameters
            class(Abstract_Diameters), intent(inout) :: this
            real(DP), intent(in) :: diameter
        end subroutine Abstract_Diameters_add
        
        subroutine Abstract_Diameters_remove(this, i_particle)
        import :: DP, Abstract_Diameters
            class(Abstract_Diameters), intent(inout) :: this
            integer, intent(in) :: i_particle
        end subroutine Abstract_Diameters_remove
        
    end interface

    type, extends(Abstract_Diameters), public :: Uniform_Diameters
    private
        real(DP) :: diameter
    contains
        procedure :: set => Uniform_Diameters_set
        procedure :: get => Uniform_Diameters_get
        procedure :: add => Uniform_Diameters_add
        procedure :: remove => Uniform_Diameters_remove
    end type Uniform_Diameters

contains

!implementation Abstract_Diameters

    subroutine Abstract_Diameters_construct(this, particles_num)
        class(Abstract_Diameters), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_num
        
        this%particles_num%ptr => particles_num
    end subroutine Abstract_Diameters_construct
    
    subroutine Abstract_Diameters_destroy(this)
        class(Abstract_Diameters), intent(inout) :: this
        
        this%particles_num%ptr => null()
    end subroutine Abstract_Diameters_destroy

!end implementation Abstract_Diameters

!implementation Uniform_Diameters

    subroutine Uniform_Diameters_set(this, i_particle, diameter)
        class(Uniform_Diameters), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: diameter
        
        if (i_particle < 1 .or. this%particles_num%ptr%get() < i_particle) then
            call error_exit("Uniform_Diameters: i_particle is out of range.")
        end if
        if (diameter < 0._DP) then
            call error_exit("Uniform_Diameters: diameter is negative.")
        end if
        if (diameter < real_zero) then
            call warning_continue("Uniform_Diameters: diameter may be too small.")
        end if
        
        this%diameter = diameter
    end subroutine Uniform_Diameters_set

    pure function Uniform_Diameters_get(this, i_particle) result(diameter)
        class(Uniform_Diameters), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: diameter
        
        diameter = this%diameter
    end function Uniform_Diameters_get
    
    subroutine Uniform_Diameters_add(this, diameter)
         class(Uniform_Diameters), intent(inout) :: this
         real(DP), intent(in) :: diameter
         
    end subroutine Uniform_Diameters_add
    
    subroutine Uniform_Diameters_remove(this, i_particle)
         class(Uniform_Diameters), intent(inout) :: this
         integer, intent(in) :: i_particle
         
    end subroutine Uniform_Diameters_remove
    
!end implementation Uniform_Diameters

end module class_diameters
