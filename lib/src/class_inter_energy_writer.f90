module class_inter_energy_writer

use data_constants, only: max_line_length
use procedures_checks, only: check_string_not_empty
use class_number_to_string, only: Abstract_Number_to_String, &
    Concrete_Number_to_String, Null_Number_to_String
use module_component_energy, only: Concrete_Inter_Energy

implicit none

private

    type, public :: Concrete_Inter_Energy_Writer_Selector
        logical :: write_long
    end type Concrete_Inter_Energy_Writer_Selector

    type, abstract, public :: Abstract_Inter_Energy_Writer
    private
        integer :: unit
        type(Concrete_Number_to_String) :: string_short
        class(Abstract_Number_to_String), allocatable :: string_long
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

    subroutine Abstract_Inter_Energy_Writer_construct(this, filename, energy_selector)
        class(Abstract_Inter_Energy_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        type(Concrete_Inter_Energy_Writer_Selector), intent(in) :: energy_selector

        character(len=:), allocatable :: legend

        call check_string_not_empty("Abstract_Inter_Energy_Writer_construct: filename", &
            filename)
        open(newunit=this%unit, recl=max_line_length, file=filename, action="write")
        legend = "# i_step    short"
        if (energy_selector%write_long) then
            allocate(Concrete_Number_to_String :: this%string_long)
            legend = legend//"    long"
        else
            allocate(Null_Number_to_String :: this%string_long)
        end if
        write(this%unit, *) legend
        deallocate(legend)
    end subroutine Abstract_Inter_Energy_Writer_construct

    subroutine Abstract_Inter_Energy_Writer_destroy(this)
        class(Abstract_Inter_Energy_Writer), intent(inout) :: this

        if (allocated(this%string_long)) deallocate(this%string_long)
        close(this%unit)
    end subroutine Abstract_Inter_Energy_Writer_destroy

    subroutine Abstract_Inter_Energy_Writer_write(this, i_step, energy)
        class(Abstract_Inter_Energy_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Concrete_Inter_Energy), intent(in) :: energy

        write(this%unit, *) i_step, this%string_short%get(energy%short), &
            this%string_long%get(energy%long)
    end subroutine Abstract_Inter_Energy_Writer_write

!end implementation Abstract_Inter_Energy_Writer

!implementation Null_Inter_Energy_Writer

    subroutine Null_Inter_Energy_Writer_construct(this, filename, energy_selector)
        class(Null_Inter_Energy_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        type(Concrete_Inter_Energy_Writer_Selector), intent(in) :: energy_selector
    end subroutine Null_Inter_Energy_Writer_construct

    subroutine Null_Inter_Energy_Writer_destroy(this)
        class(Null_Inter_Energy_Writer), intent(inout) :: this
    end subroutine Null_Inter_Energy_Writer_destroy

    subroutine Null_Inter_Energy_Writer_write(this, i_step, energy)
        class(Null_Inter_Energy_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Concrete_Inter_Energy), intent(in) :: energy
    end subroutine Null_Inter_Energy_Writer_write

!end implementation Null_Inter_Energy_Writer

end module class_inter_energy_writer
