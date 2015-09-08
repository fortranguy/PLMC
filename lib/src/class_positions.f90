module class_positions

use, intrinsic :: iso_fortran_env, only: DP => REAL64, error_unit
use data_geometry, only: num_dimensions
use data_precisions, only: real_zero
use module_error, only: error_exit
use class_particles_number, only: Abstract_Particles_Number, Abstract_Particles_Number_Pointer
use procedures_coordinates, only: increase_coordinates_size

implicit none

private

    type, abstract, public :: Abstract_Positions
    private
        type(Abstract_Particles_Number_Pointer) :: particles_num
        real(DP), allocatable :: positions(:, :)
    contains
        procedure :: construct => Abstract_Positions_construct
        procedure :: destroy => Abstract_Positions_destroy
            
        procedure :: set => Abstract_Positions_set
        procedure :: get => Abstract_Positions_get
        procedure :: add => Abstract_Positions_add
        procedure :: remove => Abstract_Positions_remove
    end type Abstract_Positions
    
    type, extends(Abstract_Positions), public :: Concrete_Positions
    
    end type Concrete_Positions
    
contains

!implementation Abstract_Positions

    subroutine Abstract_Positions_construct(this, particles_num)
        class(Abstract_Positions), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_num
        
        this%particles_num%ptr => particles_num
        if (this%particles_num%ptr%get() == 0) then
            allocate(this%positions(num_dimensions, 1))
        else
            allocate(this%positions(num_dimensions, this%particles_num%ptr%get()))
        end if
        this%positions = 0._DP
    end subroutine Abstract_Positions_construct
        
    subroutine Abstract_Positions_destroy(this)
        class(Abstract_Positions), intent(inout) :: this
        
        if (allocated(this%positions)) deallocate(this%positions)
        this%particles_num%ptr => null()
    end subroutine Abstract_Positions_destroy
    
    subroutine Abstract_Positions_set(this, i_particle, position)
        class(Abstract_Positions), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: position(:)
        
        if (size(position) /= num_dimensions) then
            call error_exit("Abstract_Positions: wrong number of dimensions (size).")
        end if
        this%positions(:, i_particle) = position
    end subroutine Abstract_Positions_set
    
    pure function Abstract_Positions_get(this, i_particle) result(position)
        class(Abstract_Positions), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: position(num_dimensions)
        
        position = this%positions(:, i_particle)
    end function Abstract_Positions_get
    
    subroutine Abstract_Positions_add(this, position)
        class(Abstract_Positions), intent(inout) :: this
        real(DP), intent(in) :: position(num_dimensions)
        
        if (size(this%positions, 2) < this%particles_num%ptr%get()) then
            call increase_coordinates_size(this%positions)
        end if
        call this%set(this%particles_num%ptr%get(), position)
    end subroutine Abstract_Positions_add
    
    subroutine Abstract_Positions_remove(this, i_particle)
        class(Abstract_Positions), intent(inout) :: this
        integer, intent(in) :: i_particle
        
        if (i_particle < 1 .or. this%particles_num%ptr%get()+1 < i_particle) then
            call error_exit("Uniform_Spheres: i_particle is out of range.")
        end if        
        if (i_particle < this%particles_num%ptr%get()+1) then
            call this%set(i_particle, this%get(this%particles_num%ptr%get()+1))
        end if
    end subroutine Abstract_Positions_remove
    
!end implementation Abstract_Positions

end module class_positions
