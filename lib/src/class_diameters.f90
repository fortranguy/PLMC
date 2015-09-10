module class_diameters

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use procedures_errors, only: warning_continue
use procedures_checks, only: check_in_range, check_positive
use class_particles_number, only: Abstract_Particles_Number

implicit none

private

    type, abstract, public :: Abstract_Diameters
        class(Abstract_Particles_Number), pointer :: particles_number
    contains
        procedure :: construct => Abstract_Diameters_construct
        procedure :: destroy => Abstract_Diameters_destroy
        procedure(Abstract_Diameters_set), deferred :: set
        procedure(Abstract_Diameters_get), deferred :: get
        procedure(Abstract_Diameters_add), deferred :: add
        procedure(Abstract_Diameters_remove), deferred :: remove
    end type Abstract_Diameters
    
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

    subroutine Abstract_Diameters_construct(this, particles_number)
        class(Abstract_Diameters), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_number
        
        this%particles_number => particles_number
    end subroutine Abstract_Diameters_construct
    
    subroutine Abstract_Diameters_destroy(this)
        class(Abstract_Diameters), intent(inout) :: this
        
        this%particles_number => null()
    end subroutine Abstract_Diameters_destroy

!end implementation Abstract_Diameters

!implementation Uniform_Diameters

    subroutine Uniform_Diameters_set(this, i_particle, diameter)
        class(Uniform_Diameters), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: diameter
        
        call check_in_range("Uniform_Diameter", this%particles_number%get(), &
                            "i_particle", i_particle)
        call check_positive("Uniform_Diameters", "diameter", diameter)
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
        
        if (abs(diameter - this%diameter) > real_zero) then
            call warning_continue("Uniform_Diameters: adding diameter is different.")
        end if
        call this%set(this%particles_number%get(), diameter)
    end subroutine Uniform_Diameters_add
    
    subroutine Uniform_Diameters_remove(this, i_particle)
        class(Uniform_Diameters), intent(inout) :: this
        integer, intent(in) :: i_particle

        call check_in_range("Uniform_Diameters", this%particles_number%get(), "i_particle", i_particle)
    end subroutine Uniform_Diameters_remove
    
!end implementation Uniform_Diameters

end module class_diameters
