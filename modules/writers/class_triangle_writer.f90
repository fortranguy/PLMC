module class_triangle_writer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: max_line_length
use procedures_checks, only: check_string_not_empty
use class_number_to_string, only: Concrete_Number_to_String, Null_Number_to_String
use types_line_observables, only: Concrete_Line_Observables
use module_string_wrapper, only: String_Wrapper, strings_wrapper_destroy

implicit none

private

    type, public :: Concrete_Line_Selector
        logical, allocatable :: line(:)
    end type Concrete_Line_Selector

    type :: Strings_Wrapper
        type(String_Wrapper), allocatable :: with_component(:)
    end type Strings_Wrapper

    type, abstract, public :: Abstract_Triangle_Writer
    private
        integer :: file_unit
        type(Strings_Wrapper), allocatable :: strings(:)
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

    subroutine Abstract_construct(this, filename, selector)
        class(Abstract_Triangle_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        type(Concrete_Line_Selector), intent(in) :: selector(:)

        character(len=:), allocatable :: legend
        integer :: file_unit !strange gfortran behaviour: otherwise writes to output_unit.

        call check_string_not_empty("Abstract_Triangle_Writer: construct: filename", &
            filename)
        open(newunit=file_unit, recl=max_line_length, file=filename, action="write")
        this%file_unit = file_unit
        legend = "# i_step"
        call this%allocate_string(legend, selector)
        write(this%file_unit, *) legend
    end subroutine Abstract_construct

    subroutine Abstract_allocate_strings(this, legend, selector)
        class(Abstract_Triangle_Writer), intent(inout) :: this
        character(len=:), allocatable, intent(inout) :: legend
        type(Concrete_Line_Selector), intent(in) :: selector(:)

        type(Concrete_Number_to_String) :: string
        integer :: i_component, j_component

        allocate(this%strings(size(selector)))
        do j_component = 1, size(this%strings)
            allocate(this%strings(j_component)%with_component(j_component))
            do i_component = 1, size(this%strings(j_component)%with_component)
                if (selector(j_component)%line(i_component)) then
                    allocate(Concrete_Number_to_String :: this%strings(j_component)%&
                        with_component(i_component)%string)
                    legend = legend//"    "//string%get(i_component)//"<->"//string%get(j_component)
                else
                    allocate(Null_Number_to_String :: this%strings(j_component)%&
                        with_component(i_component)%string)
                end if
            end do
        end do
    end subroutine Abstract_allocate_strings

    subroutine Abstract_destroy(this)
        class(Abstract_Triangle_Writer), intent(inout) :: this

        integer :: i_component

        if (allocated(this%strings)) then
            do i_component = size(this%strings), 1, -1
                call strings_wrapper_destroy(this%strings(i_component)%with_component)
            end do
            deallocate(this%strings)
        end if
        close(this%file_unit)
    end subroutine Abstract_destroy

    subroutine Abstract_write(this, i_step, triangle)
        class(Abstract_Triangle_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Concrete_Line_Observables), intent(in) :: triangle(:)

        character(len=:), allocatable :: string
        integer :: i_component, j_component

        string = ""
        do j_component = 1, size(this%strings)
            do i_component = 1, size(this%strings(j_component)%with_component)
                associate(string_ij => this%strings(j_component)%with_component(i_component)%&
                    string, energy_ij => triangle(j_component)%&
                    line(i_component))
                    string = string//string_ij%get(energy_ij)
                end associate
            end do
        end do
        write(this%file_unit, *) i_step, string
    end subroutine Abstract_write

!end implementation Abstract_Triangle_Writer

!implementation Null_Triangle_Writer

    subroutine Null_construct(this, filename, selector)
        class(Null_Triangle_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        type(Concrete_Line_Selector), intent(in) :: selector(:)
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Triangle_Writer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_write(this, i_step, triangle)
        class(Null_Triangle_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Concrete_Line_Observables), intent(in) :: triangle(:)
    end subroutine Null_write

!end implementation Null_Triangle_Writer

end module class_triangle_writer
