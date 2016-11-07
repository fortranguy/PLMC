module classes_box_size_memento

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use classes_periodic_box, only: Abstract_Periodic_Box

implicit none

private

    type, abstract, public :: Abstract_Box_Size_Memento
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        real(DP) :: box_size(num_dimensions) = 0._DP
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: target => Abstract_target
        procedure :: save => Abstract_save
        procedure :: get => Abstract_get
    end type Abstract_Box_Size_Memento

    type, extends(Abstract_Box_Size_Memento), public :: Retentive_Box_Size_Memento

    end type Retentive_Box_Size_Memento

    type, extends(Abstract_Box_Size_Memento), public :: Forgetful_Box_Size_Memento
    contains
        procedure :: construct => Forgetful_construct
        procedure :: destroy => Forgetful_destroy
        procedure :: save => Forgetful_save
        procedure :: get => Forgetful_get
    end type Forgetful_Box_Size_Memento

    type, extends(Abstract_Box_Size_Memento), public :: Null_Box_Size_Memento
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: target => Null_target
        procedure :: save => Null_save
        procedure :: get => Null_get
    end type Null_Box_Size_Memento

contains

!implementation Abstract_Box_Size_Memento

    subroutine Abstract_construct(this, periodic_box)
        class(Abstract_Box_Size_Memento), intent(out) :: this
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        call this%target(periodic_box)
        call this%save()
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Box_Size_Memento), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_target(this, periodic_box)
        class(Abstract_Box_Size_Memento), intent(inout) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box

        this%periodic_box => periodic_box
    end subroutine Abstract_target

    subroutine Abstract_save(this)
        class(Abstract_Box_Size_Memento), intent(inout) :: this

        this%box_size = this%periodic_box%get_size()
    end subroutine Abstract_save

    pure function Abstract_get(this) result(box_size)
        real(DP) :: box_size(num_dimensions)
        class(Abstract_Box_Size_Memento), intent(in) :: this

        box_size = this%box_size
    end function Abstract_get

!end implementation Abstract_Box_Size_Memento

!implementation Forgetful_Box_Size_Memento

    subroutine Forgetful_construct(this, periodic_box)
        class(Forgetful_Box_Size_Memento), intent(out) :: this
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        call this%target(periodic_box)
    end subroutine Forgetful_construct

    subroutine Forgetful_destroy(this)
        class(Forgetful_Box_Size_Memento), intent(inout) :: this
    end subroutine Forgetful_destroy

    subroutine Forgetful_save(this)
        class(Forgetful_Box_Size_Memento), intent(inout) :: this
    end subroutine Forgetful_save

    pure function Forgetful_get(this) result(box_size)
        real(DP) :: box_size(num_dimensions)
        class(Forgetful_Box_Size_Memento), intent(in) :: this

        box_size = this%periodic_box%get_size()
    end function Forgetful_get

!end implementation Forgetful_Box_Size_Memento

!implementation Null_Box_Size_Memento

    subroutine Null_construct(this, periodic_box)
        class(Null_Box_Size_Memento), intent(out) :: this
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Box_Size_Memento), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_target(this, periodic_box)
        class(Null_Box_Size_Memento), intent(inout) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
    end subroutine Null_target

    subroutine Null_save(this)
        class(Null_Box_Size_Memento), intent(inout) :: this
    end subroutine Null_save

    pure function Null_get(this) result(box_size)
        real(DP) :: box_size(num_dimensions)
        class(Null_Box_Size_Memento), intent(in) :: this
        box_size = 0._DP
    end function Null_get

!end implementation Null_Box_Size_Memento

end module classes_box_size_memento
