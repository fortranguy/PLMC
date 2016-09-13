module classes_box_volume_memento

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_periodic_box, only: Abstract_Periodic_Box

implicit none

private

    type, abstract, public :: Abstract_Box_Volume_Memento
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        real(DP) :: box_volume = 0._DP
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: save => Abstract_save
        procedure :: get_ratio => Abstract_get_ratio
    end type Abstract_Box_Volume_Memento

    type, extends(Abstract_Box_Volume_Memento), public :: Concrete_Box_Volume_Memento

    end type Concrete_Box_Volume_Memento

    type, extends(Abstract_Box_Volume_Memento), public :: Null_Box_Volume_Memento
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: save => Null_save
        procedure :: get_ratio => Null_get_ratio
    end type Null_Box_Volume_Memento

contains

!implementation Abstract_Box_Volume_Memento

    subroutine Abstract_construct(this, periodic_box)
        class(Abstract_Box_Volume_Memento), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box

        this%periodic_box => periodic_box
        call this%save()
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Box_Volume_Memento), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_save(this)
        class(Abstract_Box_Volume_Memento), intent(inout) :: this

        this%box_volume = product(this%periodic_box%get_size())
    end subroutine Abstract_save

    !> \[ \frac{V_\text{saved}}{V} \]
    pure real(DP) function Abstract_get_ratio(this) result(ratio)
        class(Abstract_Box_Volume_Memento), intent(in) :: this

        ratio = this%box_volume / product(this%periodic_box%get_size())
    end function Abstract_get_ratio

!end implementation Abstract_Box_Volume_Memento

!implementation Null_Box_Volume_Memento

    subroutine Null_construct(this, periodic_box)
        class(Null_Box_Volume_Memento), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Box_Volume_Memento), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_save(this)
        class(Null_Box_Volume_Memento), intent(inout) :: this
    end subroutine Null_save

    pure real(DP) function Null_get_ratio(this) result(ratio)
        class(Null_Box_Volume_Memento), intent(in) :: this
        ratio = 1._DP
    end function Null_get_ratio

!end implementation Null_Box_Volume_Memento

end module classes_box_volume_memento
