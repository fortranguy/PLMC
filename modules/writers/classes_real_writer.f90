module classes_real_writer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_strings, only: max_line_length
use procedures_checks, only: check_string_not_empty

implicit none

private

    type, abstract, public :: Abstract_Real_Writer
    private
        integer :: file_unit = 0
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        generic :: write => write_scalar, write_array
        procedure, private :: write_scalar => Abstract_write_scalar
        procedure, private :: write_array => Abstract_write_array
    end type Abstract_Real_Writer

    type, extends(Abstract_Real_Writer), public :: Concrete_Real_Writer

    end type Concrete_Real_Writer

    type, extends(Abstract_Real_Writer), public :: Null_Real_Writer
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: write_scalar => Null_write_scalar
        procedure :: write_array => Null_write_array
    end type Null_Real_Writer

contains

!implementation Abstract_Real_Writer

    subroutine Abstract_construct(this, filename)
        class(Abstract_Real_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename

        character(len=:), allocatable :: legend
        integer :: file_unit

        call check_string_not_empty("Abstract_Real_Writer: construct: filename", filename)
        open(newunit=file_unit, recl=max_line_length, file=filename, action="write")
        this%file_unit = file_unit
        legend = "# i_step  observable"
        write(this%file_unit, *) legend
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Real_Writer), intent(inout) :: this

        close(this%file_unit)
    end subroutine Abstract_destroy

    subroutine Abstract_write_scalar(this, i_step, observable)
        class(Abstract_Real_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: observable

        write(this%file_unit, *) i_step, observable
    end subroutine Abstract_write_scalar

    !> This subroutine differs from [[classes_line_writer:Abstract_write_reals]]:
    !> there is no local selection (filter).
    !> @todo Merge them?
    subroutine Abstract_write_array(this, i_step, observable)
        class(Abstract_Real_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: observable(:)

        write(this%file_unit, *) i_step, observable
    end subroutine Abstract_write_array

!end implementation Abstract_Real_Writer

!implementation Null_Real_Writer

    subroutine Null_construct(this, filename)
        class(Null_Real_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Real_Writer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_write_scalar(this, i_step, observable)
        class(Null_Real_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: observable
    end subroutine Null_write_scalar

    subroutine Null_write_array(this, i_step, observable)
        class(Null_Real_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: observable(:)
    end subroutine Null_write_array

!end implementation Null_Real_Writer

end module classes_real_writer
