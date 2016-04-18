module class_energy_writer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_strings, only: max_line_length
use procedures_checks, only: check_string_not_empty

implicit none

private

    type, abstract, public :: Abstract_Energy_Writer
    private
        integer :: file_unit = 0
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: write => Abstract_write
    end type Abstract_Energy_Writer

    type, extends(Abstract_Energy_Writer), public :: Concrete_Energy_Writer

    end type Concrete_Energy_Writer

    type, extends(Abstract_Energy_Writer), public :: Null_Energy_Writer
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: write => Null_write
    end type Null_Energy_Writer

contains

!implementation Abstract_Energy_Writer

    subroutine Abstract_construct(this, filename)
        class(Abstract_Energy_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename

        character(len=:), allocatable :: legend
        integer :: file_unit

        call check_string_not_empty("Abstract_Energy_Writer: construct: filename", filename)
        open(newunit=file_unit, recl=max_line_length, file=filename, action="write")
        this%file_unit = file_unit
        legend = "# i_step  energy"
        write(this%file_unit, *) legend
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Energy_Writer), intent(inout) :: this

        close(this%file_unit)
    end subroutine Abstract_destroy

    subroutine Abstract_write(this, i_step, energy)
        class(Abstract_Energy_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: energy

        write(this%file_unit, *) i_step, energy
    end subroutine Abstract_write

!end implementation Abstract_Energy_Writer

!implementation Null_Energy_Writer

    subroutine Null_construct(this, filename)
        class(Null_Energy_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Energy_Writer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_write(this, i_step, energy)
        class(Null_Energy_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        real(DP), intent(in) :: energy
    end subroutine Null_write

!end implementation Null_Energy_Writer

end module class_energy_writer
