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
        procedure :: get => Abstract_get
    end type Abstract_Box_Volume_Memento

    type, extends(Abstract_Box_Volume_Memento), public :: Concrete_Box_Volume_Memento

    end type Concrete_Box_Volume_Memento

    type, extends(Abstract_Box_Volume_Memento), public :: Forgetful_Box_Volume_Memento
    contains
        procedure :: construct => Forgetful_construct
        procedure :: destroy => Forgetful_destroy
        procedure :: save => Forgetful_save
        procedure :: get => Forgetful_get
    end type Forgetful_Box_Volume_Memento

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

    !> \[ V_\text{s} \]
    pure real(DP) function Abstract_get(this) result(box_volume)
        class(Abstract_Box_Volume_Memento), intent(in) :: this

        box_volume = this%box_volume
    end function Abstract_get

!end implementation Abstract_Box_Volume_Memento

!implementation Forgetful_Box_Volume_Memento

    subroutine Forgetful_construct(this, periodic_box)
        class(Forgetful_Box_Volume_Memento), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
    end subroutine Forgetful_construct

    subroutine Forgetful_destroy(this)
        class(Forgetful_Box_Volume_Memento), intent(inout) :: this
    end subroutine Forgetful_destroy

    subroutine Forgetful_save(this)
        class(Forgetful_Box_Volume_Memento), intent(inout) :: this
    end subroutine Forgetful_save

    !> \[ V \]
    pure real(DP) function Forgetful_get(this) result(box_volume)
        class(Forgetful_Box_Volume_Memento), intent(in) :: this
        box_volume = product(this%periodic_box%get_size())
    end function Forgetful_get

!end implementation Forgetful_Box_Volume_Memento

end module classes_box_volume_memento
