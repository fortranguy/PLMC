module classes_box_size_writer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_periodic_box, only: Abstract_Periodic_Box

implicit none

private

    type, abstract, public :: Abstract_Box_Size_Writer
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
    contains
        procedure :: construct_new => Abstract_construct_new
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: write_new => Abstract_write_new
        procedure :: write => Abstract_write
    end type Abstract_Box_Size_Writer

    type, extends(Abstract_Box_Size_Writer), public :: Concrete_Box_Size_Writer

    end type Concrete_Box_Size_Writer

    type, extends(Abstract_Box_Size_Writer), public :: Null_Box_Size_Writer
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: write => Null_write
    end type Null_Box_Size_Writer

contains

!implementation Abstract_Box_Size_Writer

    subroutine Abstract_construct_new(this, periodic_box)
        class(Abstract_Box_Size_Writer), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box

        this%periodic_box => periodic_box
    end subroutine Abstract_construct_new

    subroutine Abstract_construct(this, periodic_box, period, basename)
        class(Abstract_Box_Size_Writer), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        integer, intent(in) :: period
        character(len=*), intent(in) :: basename
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Box_Size_Writer), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_write_new(this, coordinates_unit)
        class(Abstract_Box_Size_Writer), intent(in) :: this
        integer, intent(in) :: coordinates_unit

        write(coordinates_unit, *) "# box_size:", this%periodic_box%get_size()
    end subroutine Abstract_write_new

    subroutine Abstract_write(this, i_step)
        class(Abstract_Box_Size_Writer), intent(in) :: this
        integer, intent(in) :: i_step
    end subroutine Abstract_write

!end implementation Abstract_Box_Size_Writer

!implementation Null_Box_Size_Writer

    subroutine Null_construct(this, periodic_box, period, basename)
        class(Null_Box_Size_Writer), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        integer, intent(in) :: period
        character(len=*), intent(in) :: basename
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Box_Size_Writer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_write(this, i_step)
        class(Null_Box_Size_Writer), intent(in) :: this
        integer, intent(in) :: i_step
    end subroutine Null_write

!end implementation Null_Box_Size_Writer

end module classes_box_size_writer
