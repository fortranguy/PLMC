module classes_rectangle_writer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_strings, only: max_line_length
use procedures_checks, only: check_string_not_empty
use classes_number_to_string, only: Concrete_Number_to_String, Null_Number_to_String, &
    Number_to_String_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Rectangle_Writer
    private
        type(Number_to_String_Wrapper), allocatable :: strings(:, :)
        integer :: file_unit = 0
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure, private :: allocate_strings => Abstract_allocate_strings
        procedure :: write => Abstract_write
    end type Abstract_Rectangle_Writer

    type, extends(Abstract_Rectangle_Writer), public :: Concrete_Rectangle_Writer

    end type Concrete_Rectangle_Writer

    type, extends(Abstract_Rectangle_Writer), public :: Null_Rectangle_Writer
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: write => Null_write
    end type Null_Rectangle_Writer

contains

!implementation Abstract_Rectangle_Writer

    subroutine Abstract_construct(this, filename, selectors)
        class(Abstract_Rectangle_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        logical, intent(in) :: selectors(:, :)

        character(len=:), allocatable :: legend
        integer :: file_unit

        call check_string_not_empty("Abstract_Rectangle_Writer: construct: filename", filename)
        open(newunit=file_unit, recl=max_line_length, file=filename, action="write")
        this%file_unit = file_unit
        legend = "# i_step"
        call this%allocate_strings(legend, selectors)
        write(this%file_unit, *) legend
    end subroutine Abstract_construct

    subroutine Abstract_allocate_strings(this, legend, selectors)
        class(Abstract_Rectangle_Writer), intent(inout) :: this
        character(len=:), allocatable, intent(inout) :: legend
        logical, intent(in) :: selectors(:, :)

        type(Concrete_Number_to_String) :: string
        integer :: i_element, j_element

        allocate(this%strings(size(selectors, 1), size(selectors, 2)))
        do j_element = 1, size(this%strings, 2)
            do i_element = 1, size(this%strings, 1)
                if (selectors(i_element, j_element)) then
                    allocate(Concrete_Number_to_String :: this%strings(i_element, j_element)%&
                        string)
                    legend = legend//"    "//string%get(i_element)//"->"//string%get(j_element)
                else
                    allocate(Null_Number_to_String :: this%strings(i_element, j_element)%string)
                end if
            end do
        end do
    end subroutine Abstract_allocate_strings

    subroutine Abstract_destroy(this)
        class(Abstract_Rectangle_Writer), intent(inout) :: this

        if (allocated(this%strings)) deallocate(this%strings)
        close(this%file_unit)
    end subroutine Abstract_destroy

    subroutine Abstract_write(this, i_step, observables)
        class(Abstract_Rectangle_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: observables(:, :)

        character(len=:), allocatable :: string
        integer :: i_element, j_element

        string = ""
        do j_element = 1, size(this%strings, 2)
            do i_element = 1, size(this%strings, 1)
                string = string//this%strings(i_element, j_element)%string%&
                    get(observables(i_element, j_element))
            end do
        end do
        write(this%file_unit, *) i_step, string
    end subroutine Abstract_write

!end implementation Abstract_Rectangle_Writer

!implementation Null_Rectangle_Writer

    subroutine Null_construct(this, filename, selectors)
        class(Null_Rectangle_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        logical, intent(in) :: selectors(:, :)
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Rectangle_Writer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_write(this, i_step, observables)
        class(Null_Rectangle_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: observables(:, :)
    end subroutine Null_write

!end implementation Null_Rectangle_Writer

end module classes_rectangle_writer
