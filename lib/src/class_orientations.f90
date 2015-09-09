module class_orientations

use, intrinsic :: iso_fortran_env, only: DP => REAL64, error_unit
use data_geometry, only: num_dimensions
use data_precisions, only: real_zero
use procedures_errors, only: error_exit, warning_continue
use class_particles_number, only: Abstract_Particles_Number
use procedures_coordinates, only: increase_coordinates_size

implicit none

private

    type, abstract, public :: Abstract_Orientations
    contains
        procedure(Abstract_Orientations_construct), deferred :: construct
        procedure(Abstract_Orientations_destroy), deferred :: destroy
            
        procedure(Abstract_Orientations_set), deferred :: set
        procedure(Abstract_Orientations_get), deferred :: get
        procedure(Abstract_Orientations_add), deferred :: add
        procedure(Abstract_Orientations_remove), deferred :: remove
    end type Abstract_Orientations
    
    abstract interface
    
        subroutine Abstract_Orientations_construct(this, particles_num)
        import :: Abstract_Particles_Number, Abstract_Orientations
            class(Abstract_Orientations), intent(out) :: this
            class(Abstract_Particles_Number), target, intent(in) :: particles_num
        end subroutine Abstract_Orientations_construct
        
        subroutine Abstract_Orientations_destroy(this)
        import :: Abstract_Orientations
            class(Abstract_Orientations), intent(inout) :: this
        end subroutine Abstract_Orientations_destroy
        
        subroutine Abstract_Orientations_set(this, i_particle, orientation)
        import :: DP, Abstract_Orientations
            class(Abstract_Orientations), intent(inout) :: this
            integer, intent(in) :: i_particle
            real(DP), intent(in) :: orientation(:)
            
        end subroutine Abstract_Orientations_set
        
        pure function Abstract_Orientations_get(this, i_particle) result(orientation)
        import :: DP, num_dimensions, Abstract_Orientations
            class(Abstract_Orientations), intent(in) :: this
            integer, intent(in) :: i_particle
            real(DP) :: orientation(num_dimensions)
            
        end function Abstract_Orientations_get
        
        subroutine Abstract_Orientations_add(this, orientation)
        import :: DP, num_dimensions, Abstract_Orientations
            class(Abstract_Orientations), intent(inout) :: this
            real(DP), intent(in) :: orientation(num_dimensions)
            
        end subroutine Abstract_Orientations_add
        
        subroutine Abstract_Orientations_remove(this, i_particle)
        import :: Abstract_Orientations
            class(Abstract_Orientations), intent(inout) :: this
            integer, intent(in) :: i_particle
            
        end subroutine Abstract_Orientations_remove
        
    end interface
    
    type, extends(Abstract_Orientations), public :: Null_Orientations
    contains
        procedure :: construct => Null_Orientations_construct
        procedure :: destroy => Null_Orientations_destroy
            
        procedure :: set => Null_Orientations_set
        procedure :: get => Null_Orientations_get
        procedure :: add => Null_Orientations_add
        procedure :: remove => Null_Orientations_remove
    end type
    
    type, extends(Abstract_Orientations), public :: Concrete_Orientations
    private
        class(Abstract_Particles_Number), pointer :: particles_num
        real(DP), allocatable :: orientations(:, :)
    contains
        procedure :: construct => Concrete_Orientations_construct
        procedure :: destroy => Concrete_Orientations_destroy
            
        procedure :: set => Concrete_Orientations_set
        procedure :: get => Concrete_Orientations_get
        procedure :: add => Concrete_Orientations_add
        procedure :: remove => Concrete_Orientations_remove
    end type Concrete_Orientations
    
contains

!implementation Null_Orientations

    subroutine Null_Orientations_construct(this, particles_num)
        class(Null_Orientations), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_num
        
    end subroutine Null_Orientations_construct
    
    subroutine Null_Orientations_destroy(this)
        class(Null_Orientations), intent(inout) :: this
        
    end subroutine Null_Orientations_destroy
    
    subroutine Null_Orientations_set(this, i_particle, orientation)
        class(Null_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: orientation(:)
        
    end subroutine Null_Orientations_set
    
    pure function Null_Orientations_get(this, i_particle) result(orientation)
        class(Null_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: orientation(num_dimensions)
        
        orientation = 0._DP
    end function Null_Orientations_get
    
    subroutine Null_Orientations_add(this, orientation)
        class(Null_Orientations), intent(inout) :: this
        real(DP), intent(in) :: orientation(num_dimensions)
        
    end subroutine Null_Orientations_add
    
    subroutine Null_Orientations_remove(this, i_particle)
        class(Null_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        
    end subroutine Null_Orientations_remove
        
!end implementation Null_Orientations

!implementation Concrete_Orientations

    subroutine Concrete_Orientations_construct(this, particles_num)
        class(Concrete_Orientations), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_num
        
        integer :: i_particle
        
        this%particles_num => particles_num
        if (this%particles_num%get() == 0) then
            allocate(this%orientations(num_dimensions, 1))
        else
            allocate(this%orientations(num_dimensions, this%particles_num%get()))
        end if
        
        do i_particle = 1, this%particles_num%get()
            this%orientations(:, i_particle) = [1.0_DP, 0._DP, 0._DP]
        end do
    end subroutine Concrete_Orientations_construct
        
    subroutine Concrete_Orientations_destroy(this)
        class(Concrete_Orientations), intent(inout) :: this
        
        if (allocated(this%orientations)) deallocate(this%orientations)
        this%particles_num => null()
    end subroutine Concrete_Orientations_destroy
    
    subroutine Concrete_Orientations_set(this, i_particle, orientation)
        class(Concrete_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: orientation(:)
        
        if (size(orientation) /= num_dimensions) then
            call error_exit("Concrete_Orientations: wrong number of dimensions (size).")
        end if
        if (abs(norm2(orientation)-1.0_DP) > real_zero) then
            call warning_continue("Concrete_Orientations: orientation may not be normed.")
        end if
        this%orientations(:, i_particle) = orientation
    end subroutine Concrete_Orientations_set
    
    pure function Concrete_Orientations_get(this, i_particle) result(orientation)
        class(Concrete_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: orientation(num_dimensions)
        
        orientation = this%orientations(:, i_particle)
    end function Concrete_Orientations_get
    
    subroutine Concrete_Orientations_add(this, orientation)
        class(Concrete_Orientations), intent(inout) :: this
        real(DP), intent(in) :: orientation(num_dimensions)
        
        if (size(this%orientations, 2) < this%particles_num%get()) then
            call increase_coordinates_size(this%orientations)
        end if
        if (abs(norm2(orientation)-1.0_DP) > real_zero) then
            call warning_continue("Concrete_Orientations: orientation may not be normed.")
        end if
        call this%set(this%particles_num%get(), orientation)
    end subroutine Concrete_Orientations_add
    
    subroutine Concrete_Orientations_remove(this, i_particle)
        class(Concrete_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        
        if (i_particle < 1 .or. this%particles_num%get() < i_particle) then
            call error_exit("Uniform_Spheres: i_particle is out of range.")
        end if
        if (i_particle < this%particles_num%get()) then
            call this%set(i_particle, this%get(this%particles_num%get()))
        end if
    end subroutine Concrete_Orientations_remove
    
!end implementation Concrete_Orientations

end module class_orientations
