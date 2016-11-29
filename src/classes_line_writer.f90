module classes_line_writer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_strings, only: max_line_length
use procedures_checks, only: check_string_not_empty
use classes_number_to_string, only: Concrete_Number_to_String, Null_Number_to_String, &
    Number_to_String_Wrapper
use procedures_string_factory, only: string_destroy => destroy

implicit none

private

    type, abstract, public :: Abstract_Line_Writer
    private
        type(Number_to_String_Wrapper), allocatable :: strings(:)
        integer :: file_unit = 0
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        generic :: write => write_reals, write_integers
        procedure, private :: write_reals => Abstract_write_reals
        procedure, private :: write_integers => Abstract_write_integers
        procedure, private :: allocate_strings => Abstract_allocate_strings
    end type Abstract_Line_Writer

    type, extends(Abstract_Line_Writer), public :: Concrete_Line_Writer

    end type Concrete_Line_Writer

    type, extends(Abstract_Line_Writer), public :: Null_Line_Writer
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure, private :: write_reals => Null_write_reals
        procedure, private :: write_integers => Null_write_integers
    end type Null_Line_Writer

    type, public :: Line_Writer_Wrapper
        class(Abstract_Line_Writer), allocatable :: writer
    end type Line_Writer_Wrapper

contains

!implementation Abstract_Line_Writer

    subroutine Abstract_construct(this, filename, selector)
        class(Abstract_Line_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        logical, intent(in) :: selector(:)

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
        integer :: i_element

        allocate(this%strings(size(selector)))
        do i_element = 1, size(this%strings)
            if (selector(i_element)) then
                allocate(Concrete_Number_to_String :: this%strings(i_element)%string)
                legend = legend//"    "//string%get(i_element)
            else
                allocate(Null_Number_to_String :: this%strings(i_element)%string)
            end if
        end do
    end subroutine Abstract_allocate_strings

    subroutine Abstract_destroy(this)
        class(Abstract_Line_Writer), intent(inout) :: this

        call string_destroy(this%strings)
        close(this%file_unit)
    end subroutine Abstract_destroy

    subroutine Abstract_write_reals(this, i_step, observables)
        class(Abstract_Line_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: observables(:)

        character(len=:), allocatable :: string
        integer :: i_element

        string = ""
        do i_element = 1, size(this%strings)
            string = string//this%strings(i_element)%string%get(observables(i_element))
        end do
        write(this%file_unit, *) i_step, string
    end subroutine Abstract_write_reals

    subroutine Abstract_write_integers(this, i_step, observables)
        class(Abstract_Line_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        integer, intent(in) :: observables(:)

        character(len=:), allocatable :: string
        integer :: i_element

        string = ""
        do i_element = 1, size(this%strings)
            string = string//"    "//this%strings(i_element)%string%get(observables(i_element))
        end do
        write(this%file_unit, *) i_step, string
    end subroutine Abstract_write_integers

!end implementation Abstract_Line_Writer

!implementation Null_Line_Writer

    subroutine Null_construct(this, filename, selector)
        class(Null_Line_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        logical, intent(in) :: selector(:)
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Line_Writer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_write_reals(this, i_step, observables)
        class(Null_Line_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: observables(:)
    end subroutine Null_write_reals

    subroutine Null_write_integers(this, i_step, observables)
        class(Null_Line_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        integer, intent(in) :: observables(:)
    end subroutine Null_write_integers

!end implementation Null_Line_Writer

end module classes_line_writer
