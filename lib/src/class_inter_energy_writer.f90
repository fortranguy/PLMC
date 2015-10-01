module class_inter_energy_writer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use module_data, only: test_empty_string

implicit none

    type, abstract, public :: Abstract_Inter_Energy_Writer
    private
        integer :: unit
        real(DP), pointer :: energy => null()
    contains
        procedure :: construct => Abstract_Inter_Energy_Writer_construct
        procedure :: destroy => Abstract_Inter_Energy_Writer_destroy
        procedure :: write => Abstract_Inter_Energy_Writer_write
    end type Abstract_Inter_Energy_Writer

    type, extends(Abstract_Inter_Energy_Writer), public :: Concrete_Inter_Energy_Writer

    end type Concrete_Inter_Energy_Writer

    type, extends(Abstract_Inter_Energy_Writer), public :: Null_Inter_Energy_Writer
    contains
        procedure :: construct => Null_Inter_Energy_Writer_construct
        procedure :: destroy => Null_Inter_Energy_Writer_destroy
        procedure :: write => Null_Inter_Energy_Writer_write
    end type Null_Inter_Energy_Writer

contains

!implementation Abstract_Inter_Energy_Writer

    subroutine Abstract_Inter_Energy_Writer_construct(this, energy, filename)
        class(Abstract_Inter_Energy_Writer), intent(out) :: this
        real(DP), target, intent(in) :: energy
        character(len=*), intent(in) :: filename

        this%energy => energy
        call test_empty_string("Abstract_Inter_Energy_Writer_construct: filename", filename)
        open(newunit=this%unit, recl=4096, file=filename, action="write")
    end subroutine Abstract_Inter_Energy_Writer_construct

    subroutine Abstract_Inter_Energy_Writer_destroy(this)
        class(Abstract_Inter_Energy_Writer), intent(inout) :: this

        close(this%unit)
        this%energy => null()
    end subroutine Abstract_Inter_Energy_Writer_destroy

    subroutine Abstract_Inter_Energy_Writer_write(this, i_step)
        class(Abstract_Inter_Energy_Writer), intent(in) :: this
        integer, intent(in) :: i_step

        write(this%unit, *) i_step, this%energy
    end subroutine Abstract_Inter_Energy_Writer_write

!end implementation Abstract_Inter_Energy_Writer

!implementation Null_Inter_Energy_Writer

    subroutine Null_Inter_Energy_Writer_construct(this, energy, filename)
        class(Null_Inter_Energy_Writer), intent(out) :: this
        real(DP), target, intent(in) :: energy
        character(len=*), intent(in) :: filename
    end subroutine Null_Inter_Energy_Writer_construct

    subroutine Null_Inter_Energy_Writer_destroy(this)
        class(Null_Inter_Energy_Writer), intent(inout) :: this
    end subroutine Null_Inter_Energy_Writer_destroy

    subroutine Null_Inter_Energy_Writer_write(this, i_step)
        class(Null_Inter_Energy_Writer), intent(in) :: this
        integer, intent(in) :: i_step
    end subroutine Null_Inter_Energy_Writer_write

!end implementation Null_Inter_Energy_Writer

end module class_inter_energy_writer
