module classes_square_writer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_strings, only: max_line_length
use procedures_checks, only: check_string_not_empty
use classes_number_to_string, only: Concrete_Number_to_String, Null_Number_to_String
use module_string_wrapper, only: String_Wrapper, strings_wrapper_destroy

implicit none

private

    type, abstract, public :: Abstract_Square_Writer
    private
        type(String_Wrapper), allocatable :: strings(:, :)
        integer :: file_unit = 0
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure, private :: allocate_strings => Abstract_allocate_strings
        procedure :: write => Abstract_write
    end type Abstract_Square_Writer

    type, extends(Abstract_Square_Writer), public :: Concrete_Square_Writer

    end type Concrete_Square_Writer

    type, extends(Abstract_Square_Writer), public :: Null_Square_Writer
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: write => Null_write
    end type Null_Square_Writer

contains

!implementation Abstract_Square_Writer

    subroutine Abstract_construct(this, selectors, filename)
        class(Abstract_Square_Writer), intent(out) :: this
        logical, intent(in) :: selectors(:, :)
        character(len=*), intent(in) :: filename

        character(len=:), allocatable :: legend
        integer :: file_unit

        call check_string_not_empty("Abstract_Square_Writer: construct: filename", filename)
        open(newunit=file_unit, recl=max_line_length, file=filename, action="write")
        this%file_unit = file_unit
        legend = "# i_step"
        call this%allocate_strings(legend, selectors)
        write(this%file_unit, *) legend
    end subroutine Abstract_construct

    subroutine Abstract_allocate_strings(this, legend, selectors)
        class(Abstract_Square_Writer), intent(inout) :: this
        character(len=:), allocatable, intent(inout) :: legend
        logical, intent(in) :: selectors(:, :)

        type(Concrete_Number_to_String) :: string
        integer :: i_component, j_component

        allocate(this%strings(size(selectors, 1), size(selectors, 2)))
        do j_component = 1, size(this%strings, 2)
            do i_component = 1, size(this%strings, 1)
                if (selectors(i_component, j_component)) then
                    allocate(Concrete_Number_to_String :: this%strings(i_component, j_component)%&
                        string)
                    legend = legend//"    "//string%get(i_component)//"->"//string%get(j_component)
                else
                    allocate(Null_Number_to_String :: this%strings(i_component, j_component)%string)
                end if
            end do
        end do
    end subroutine Abstract_allocate_strings

    subroutine Abstract_destroy(this)
        class(Abstract_Square_Writer), intent(inout) :: this

        if (allocated(this%strings)) deallocate(this%strings)
        close(this%file_unit)
    end subroutine Abstract_destroy

    subroutine Abstract_write(this, i_step, square)
        class(Abstract_Square_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: square(:, :)

        character(len=:), allocatable :: string
        integer :: i_component, j_component

        string = ""
        do j_component = 1, size(this%strings, 2)
            do i_component = 1, size(this%strings, 1)
                string = string//this%strings(i_component, j_component)%string%&
                    get(square(i_component, j_component))
            end do
        end do
        write(this%file_unit, *) i_step, string
    end subroutine Abstract_write

!end implementation Abstract_Square_Writer

!implementation Null_Square_Writer

    subroutine Null_construct(this, selectors, filename)
        class(Null_Square_Writer), intent(out) :: this
        logical, intent(in) :: selectors(:, :)
        character(len=*), intent(in) :: filename
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Square_Writer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_write(this, i_step, square)
        class(Null_Square_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: square(:, :)
    end subroutine Null_write

!end implementation Null_Square_Writer

end module classes_square_writer
