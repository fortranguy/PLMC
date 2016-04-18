module class_line_writer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_strings, only: max_line_length
use procedures_checks, only: check_string_not_empty
use class_number_to_string, only: Concrete_Number_to_String, Null_Number_to_String
use module_string_wrapper, only: String_Wrapper, strings_wrapper_destroy

implicit none

private

    type, abstract, public :: Abstract_Line_Writer
    private
        type(String_Wrapper), allocatable :: strings(:)
        integer :: file_unit = 0
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: write => Abstract_write
        procedure, private :: allocate_strings => Abstract_allocate_strings
    end type Abstract_Line_Writer

    type, extends(Abstract_Line_Writer), public :: Concrete_Line_Writer

    end type Concrete_Line_Writer

    type, extends(Abstract_Line_Writer), public :: Null_Line_Writer
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: write => Null_write
    end type Null_Line_Writer

contains

!implementation Abstract_Line_Writer

    subroutine Abstract_construct(this, selector, filename)
        class(Abstract_Line_Writer), intent(out) :: this
        logical, intent(in) :: selector(:)
        character(len=*), intent(in) :: filename

        character(len=:), allocatable :: legend
        integer :: file_unit

        call check_string_not_empty("Abstract_Line_Writer: construct: filename", filename)
        open(newunit=file_unit, recl=max_line_length, file=filename, action="write")
        this%file_unit = file_unit
        legend = "# i_step"
        call this%allocate_strings(legend, selector)
        write(this%file_unit, *) legend
    end subroutine Abstract_construct

    subroutine Abstract_allocate_strings(this, legend, selector)
        class(Abstract_Line_Writer), intent(inout) :: this
        character(len=:), allocatable, intent(inout) :: legend
        logical, intent(in) :: selector(:)

        type(Concrete_Number_to_String) :: string
        integer :: i_component

        allocate(this%strings(size(selector)))
        do i_component = 1, size(this%strings)
            if (selector(i_component)) then
                allocate(Concrete_Number_to_String :: this%strings(i_component)%string)
                legend = legend//"    "//string%get(i_component)
            else
                allocate(Null_Number_to_String :: this%strings(i_component)%string)
            end if
        end do
    end subroutine Abstract_allocate_strings

    subroutine Abstract_destroy(this)
        class(Abstract_Line_Writer), intent(inout) :: this

        call strings_wrapper_destroy(this%strings)
        close(this%file_unit)
    end subroutine Abstract_destroy

    subroutine Abstract_write(this, i_step, line)
        class(Abstract_Line_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: line(:)

        character(len=:), allocatable :: string
        integer :: i_component

        string = ""
        do i_component = 1, size(this%strings)
            string = string//this%strings(i_component)%string%get(line(i_component))
        end do
        write(this%file_unit, *) i_step, string
    end subroutine Abstract_write

!end implementation Abstract_Line_Writer

!implementation Null_Line_Writer

    subroutine Null_construct(this, selector, filename)
        class(Null_Line_Writer), intent(out) :: this
        logical, intent(in) :: selector(:)
        character(len=*), intent(in) :: filename
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Line_Writer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_write(this, i_step, line)
        class(Null_Line_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: line(:)
    end subroutine Null_write

!end implementation Null_Line_Writer

end module class_line_writer
