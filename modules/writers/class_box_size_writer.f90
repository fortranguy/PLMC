module class_box_size_writer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: max_line_length
use procedures_checks, only: check_string_not_empty, check_positive
use class_number_to_string, only: Concrete_Number_to_String
use class_periodic_box, only: Abstract_Periodic_Box

implicit none

private

    type, abstract, public :: Abstract_Box_Size_Writer
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        character(len=:), allocatable :: basename
        character(len=:), allocatable :: legend
        integer :: period = 0
        type(Concrete_Number_to_String) :: string_step
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
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

    subroutine Abstract_construct(this, periodic_box, period, basename)
        class(Abstract_Box_Size_Writer), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        integer, intent(in) :: period
        character(len=*), intent(in) :: basename

        this%periodic_box => periodic_box
        call check_string_not_empty("Abstract_Box_Size_Writer: construct: basename", basename)
        this%basename = basename
        this%legend = "# size_x    size_y    size_z"
        call check_positive("Abstract_Box_Size_Writer: construct", "period", period)
        this%period = period
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Box_Size_Writer), intent(inout) :: this

        if (allocated(this%legend)) deallocate(this%legend)
        if (allocated(this%basename)) deallocate(this%basename)
        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_write(this, i_step)
        class(Abstract_Box_Size_Writer), intent(in) :: this
        integer, intent(in) :: i_step

        integer :: unit_i

        if (mod(i_step, this%period) == 0) then
            open(newunit=unit_i, recl=max_line_length, file=this%basename//"_"//this%string_step%&
                get(i_step)//".out", action="write")
            write(unit_i, *) this%legend
            write(unit_i, *) this%periodic_box%get_size()
            close(unit_i)
        end if
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

end module class_box_size_writer
