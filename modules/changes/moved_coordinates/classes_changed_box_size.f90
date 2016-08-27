module classes_changed_box_size

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_increase_factor
use classes_moved_coordinates, only: Abstract_Moved_Coordinates
use classes_changed_box_size_ratio, only: Abstract_Changed_Box_Size_Ratio
use module_move_tuning, only: Concrete_Move_Tuning_Parameters, set_increase_factor

implicit none

private

    type, extends(Abstract_Moved_Coordinates), abstract, public :: Abstract_Changed_Box_Size
    private
        class(Abstract_Changed_Box_Size_Ratio), allocatable :: ratio
        type(Concrete_Move_Tuning_Parameters) :: tuning_parameters
        real(DP) :: current_increase_factor = 1._DP
        logical :: max_factor_reached = .false.
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: increase_delta => Abstract_increase_delta
        procedure :: decrease_delta => Abstract_decrease_delta
        procedure :: get_ratio => Abstract_get_ratio
    end type Abstract_Changed_Box_Size

    type, extends(Abstract_Changed_Box_Size), public :: Concrete_Changed_Box_Size

    end type Concrete_Changed_Box_Size

    type, extends(Abstract_Changed_Box_Size), public :: Null_Changed_Box_Size
    contains
        procedure :: construct => Null_construct
        procedure :: increase_delta => Null_increase_delta
        procedure :: decrease_delta => Null_decrease_delta
        procedure :: get_ratio => Null_get_ratio
    end type Null_Changed_Box_Size

contains

!implementation Abstract_Changed_Box_Size

    subroutine Abstract_construct(this, ratio, tuning_parameters)
        class(Abstract_Changed_Box_Size), intent(out) :: this
        class(Abstract_Changed_Box_Size_Ratio), intent(in) :: ratio
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters

        allocate(this%ratio, source=ratio)
        call check_increase_factor("Abstract_Changed_Box_Size: construct", "increase_factor", &
            tuning_parameters%increase_factor)
        this%tuning_parameters%increase_factor = tuning_parameters%increase_factor
        this%current_increase_factor = this%tuning_parameters%increase_factor
        call check_increase_factor("Abstract_Changed_Box_Size: construct", &
            "increase_factor_max", tuning_parameters%increase_factor_max)
        this%tuning_parameters%increase_factor_max = tuning_parameters%increase_factor_max
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Changed_Box_Size), intent(inout) :: this

        if (allocated(this%ratio)) deallocate(this%ratio)
    end subroutine Abstract_destroy

    subroutine Abstract_increase_delta(this)
        class(Abstract_Changed_Box_Size), intent(inout) :: this

        if (this%max_factor_reached) return
        call set_increase_factor("Abstract_Changed_Box_Size", this%current_increase_factor, &
            this%tuning_parameters, this%max_factor_reached)
        call this%ratio%set(this%current_increase_factor * this%ratio%get_delta())
    end subroutine Abstract_increase_delta

    subroutine Abstract_decrease_delta(this)
        class(Abstract_Changed_Box_Size), intent(inout) :: this

        call this%ratio%set(this%ratio%get_delta() / this%current_increase_factor)
    end subroutine Abstract_decrease_delta

    function Abstract_get_ratio(this) result(ratio)
        class(Abstract_Changed_Box_Size), intent(in) :: this
        real(DP) :: ratio(num_dimensions)

        ratio = this%ratio%get()
    end function Abstract_get_ratio

!end implementation Abstract_Changed_Box_Size

!implementation Null_Changed_Box_Size

    subroutine Null_construct(this, ratio, tuning_parameters)
        class(Null_Changed_Box_Size), intent(out) :: this
        class(Abstract_Changed_Box_Size_Ratio), intent(in) :: ratio
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
    end subroutine Null_construct

    subroutine Null_increase_delta(this)
        class(Null_Changed_Box_Size), intent(inout) :: this
    end subroutine Null_increase_delta

    subroutine Null_decrease_delta(this)
        class(Null_Changed_Box_Size), intent(inout) :: this
    end subroutine Null_decrease_delta

    function Null_get_ratio(this) result(ratio)
        real(DP) :: ratio(num_dimensions)
        class(Null_Changed_Box_Size), intent(in) :: this
        ratio = 1._DP
    end function Null_get_ratio

!end implementation Null_Changed_Box_Size


end module classes_changed_box_size
