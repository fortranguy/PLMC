module class_small_move

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_checks, only: check_3d_array, check_positive
use class_positions, only: Abstract_Positions

implicit none

    type, abstract, public :: Abstract_Small_Move
    private
        class(Abstract_Positions), pointer :: positions
        real(DP) :: delta(num_dimensions)
    contains
        procedure :: construct => Abstract_Small_Move_construct
        procedure :: destroy => Abstract_Small_Move_destroy
        procedure :: set => Abstract_Small_Move_set
        procedure :: get => Abstract_Small_Move_get
    end type Abstract_Small_Move

    type, extends(Abstract_Small_Move), public :: Concrete_Small_Move

    end type Concrete_Small_Move
    
contains

    subroutine Abstract_Small_Move_construct(this, positions)
        class(Abstract_Small_Move), intent(out) :: this
        class(Abstract_Positions), target, intent(in) :: positions
        
        this%positions => positions
    end subroutine Abstract_Small_Move_construct

    subroutine Abstract_Small_Move_destroy(this)
        class(Abstract_Small_Move), intent(inout) :: this
        
        this%positions => null()
    end subroutine Abstract_Small_Move_destroy

    subroutine Abstract_Small_Move_set(this, delta)
        class(Abstract_Small_Move), intent(inout) :: this
        real(DP) :: delta(:)
        
        call check_3d_array("Abstract_Small_Move", "delta", delta)
        call check_positive("Abstract_Small_Move", "delta", delta)
        this%delta = delta
    end subroutine Abstract_Small_Move_set

    function Abstract_Small_Move_get(this, i_particle) result(moved_position)
        class(Abstract_Small_Move), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: moved_position(num_dimensions)

        real(DP) :: rand(num_dimensions)

        call random_number(rand)
        moved_position = this%positions%get(i_particle) + (rand - 0.5_DP) * this%delta
    end function Abstract_Small_Move_get

end module class_small_move

