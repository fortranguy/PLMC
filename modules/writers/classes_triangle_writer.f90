module classes_triangle_writer

use data_strings, only: max_line_length
use procedures_checks, only: check_string_not_empty
use classes_number_to_string, only: Concrete_Number_to_String, Null_Number_to_String
use types_reals_line, only: Reals_Line
use types_logical_line, only: Concrete_Logical_Line
use module_string_wrapper, only: Strings_Line, strings_wrapper_destroy

implicit none

private

    type, abstract, public :: Abstract_Triangle_Writer
    private
        type(Strings_Line), allocatable :: strings(:)
        integer :: file_unit = 0
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: write => Abstract_write
        procedure, private :: allocate_string => Abstract_allocate_strings
    end type Abstract_Triangle_Writer

    type, extends(Abstract_Triangle_Writer), public :: Concrete_Triangle_Writer

    end type Concrete_Triangle_Writer

    type, extends(Abstract_Triangle_Writer), public :: Null_Triangle_Writer
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: write => Null_write
    end type Null_Triangle_Writer

contains

!implementation Abstract_Triangle_Writer

    subroutine Abstract_construct(this, selectors, filename)
        class(Abstract_Triangle_Writer), intent(out) :: this
        type(Concrete_Logical_Line), intent(in) :: selectors(:)
        character(len=*), intent(in) :: filename

        character(len=:), allocatable :: legend
        integer :: file_unit !strange gfortran behaviour: otherwise writes to output_unit.

        call check_string_not_empty("Abstract_Triangle_Writer: construct: filename", filename)
        open(newunit=file_unit, recl=max_line_length, file=filename, action="write")
        this%file_unit = file_unit
        legend = "# i_step"
        call this%allocate_string(legend, selectors)
        write(this%file_unit, *) legend
    end subroutine Abstract_construct

    subroutine Abstract_allocate_strings(this, legend, selectors)
        class(Abstract_Triangle_Writer), intent(inout) :: this
        character(len=:), allocatable, intent(inout) :: legend
        type(Concrete_Logical_Line), intent(in) :: selectors(:)

        type(Concrete_Number_to_String) :: string
        integer :: i_component, j_component

        allocate(this%strings(size(selectors)))
        do j_component = 1, size(this%strings)
            allocate(this%strings(j_component)%line(j_component))
            do i_component = 1, size(this%strings(j_component)%line)
                if (selectors(j_component)%line(i_component)) then
                    allocate(Concrete_Number_to_String :: this%strings(j_component)%&
                        line(i_component)%string)
                    legend = legend//"    "//string%get(i_component)//"<->"//string%get(j_component)
                else
                    allocate(Null_Number_to_String :: this%strings(j_component)%line(i_component)%&
                        string)
                end if
            end do
        end do
    end subroutine Abstract_allocate_strings

    subroutine Abstract_destroy(this)
        class(Abstract_Triangle_Writer), intent(inout) :: this

        integer :: i_component

        if (allocated(this%strings)) then
            do i_component = size(this%strings), 1, -1
                call strings_wrapper_destroy(this%strings(i_component)%line)
            end do
            deallocate(this%strings)
        end if
        close(this%file_unit)
    end subroutine Abstract_destroy

    subroutine Abstract_write(this, i_step, observables)
        class(Abstract_Triangle_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Reals_Line), intent(in) :: observables(:)

        character(len=:), allocatable :: string
        integer :: i_component, j_component

        string = ""
        do j_component = 1, size(this%strings)
            do i_component = 1, size(this%strings(j_component)%line)
                associate(string_ij => this%strings(j_component)%line(i_component)%string, &
                    real_ij => observables(j_component)%line(i_component))
                    string = string//string_ij%get(real_ij)
                end associate
            end do
        end do
        write(this%file_unit, *) i_step, string
    end subroutine Abstract_write

!end implementation Abstract_Triangle_Writer

!implementation Null_Triangle_Writer

    subroutine Null_construct(this, selectors, filename)
        class(Null_Triangle_Writer), intent(out) :: this
        type(Concrete_Logical_Line), intent(in) :: selectors(:)
        character(len=*), intent(in) :: filename
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Triangle_Writer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_write(this, i_step, observables)
        class(Null_Triangle_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Reals_Line), intent(in) :: observables(:)
    end subroutine Null_write

!end implementation Null_Triangle_Writer

end module classes_triangle_writer
