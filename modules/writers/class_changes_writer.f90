module class_changes_writer

use data_constants, only: max_line_length
use procedures_checks, only: check_string_not_empty
use class_number_to_string, only: Abstract_Number_to_String, &
    Concrete_Number_to_String, Null_Number_to_String
use module_changes_success, only: Concrete_Changes_Success

implicit none

private

    type, public :: Concrete_Changes_Selector
        logical :: write_positions
        logical :: write_rotations
        logical :: write_exchanges
    end type Concrete_Changes_Selector

    type, abstract, public :: Abstract_Changes_Success_Writer
    private
        integer :: file_unit
        type(Concrete_Number_to_String) :: string_move
        class(Abstract_Number_to_String), allocatable :: string_rotation
        class(Abstract_Number_to_String), allocatable :: string_exchange
    contains
        procedure :: construct => Abstract_Changes_Success_Writer_construct
        procedure :: destroy => Abstract_Changes_Success_Writer_destroy
        procedure :: write => Abstract_Changes_Success_Writer_write
    end type Abstract_Changes_Success_Writer

    type, extends(Abstract_Changes_Success_Writer), public :: Concrete_Changes_Success_Writer

    end type Concrete_Changes_Success_Writer

    type, extends(Abstract_Changes_Success_Writer), public :: Null_Changes_Success_Writer
    contains
        procedure :: construct => Null_Changes_Success_Writer_construct
        procedure :: destroy => Null_Changes_Success_Writer_destroy
        procedure :: write => Null_Changes_Success_Writer_write
    end type Null_Changes_Success_Writer

contains

!implementation Abstract_Changes_Success_Writer

    subroutine Abstract_Changes_Success_Writer_construct(this, filename, changes_selector)
        class(Abstract_Changes_Success_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        type(Concrete_Changes_Selector), intent(in) :: changes_selector

        character(len=:), allocatable :: legend

        call check_string_not_empty("Abstract_Changes_Success_Writer_construct: filename", filename)
        open(newunit=this%file_unit, recl=max_line_length, file=filename, action="write")
        legend = "# i_step    moves"
        if (changes_selector%write_rotations) then
            allocate(Concrete_Number_to_String :: this%string_rotation)
            legend = legend//"    rotations"
        else
            allocate(Null_Number_to_String :: this%string_rotation)
        end if
        if (changes_selector%write_exchanges) then
            allocate(Concrete_Number_to_String :: this%string_exchange)
            legend = legend//"    exchanges"
        else
            allocate(Null_Number_to_String :: this%string_exchange)
        end if
        write(this%file_unit, *) legend
    end subroutine Abstract_Changes_Success_Writer_construct

    subroutine Abstract_Changes_Success_Writer_destroy(this)
        class(Abstract_Changes_Success_Writer), intent(inout) :: this

        if (allocated(this%string_exchange)) deallocate(this%string_exchange)
        if (allocated(this%string_rotation)) deallocate(this%string_rotation)
        close(this%file_unit)
    end subroutine Abstract_Changes_Success_Writer_destroy

    subroutine Abstract_Changes_Success_Writer_write(this, i_step, changes_success)
        class(Abstract_Changes_Success_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Concrete_Changes_Success), intent(in) :: changes_success

        write(this%file_unit, *) i_step, this%string_move%get(changes_success%move)//&
            this%string_rotation%get(changes_success%rotation)
    end subroutine Abstract_Changes_Success_Writer_write

!end implementation Abstract_Changes_Success_Writer

!implementation Null_Changes_Success_Writer

    subroutine Null_Changes_Success_Writer_construct(this, filename, changes_selector)
        class(Null_Changes_Success_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        type(Concrete_Changes_Selector), intent(in) :: changes_selector
    end subroutine Null_Changes_Success_Writer_construct

    subroutine Null_Changes_Success_Writer_destroy(this)
        class(Null_Changes_Success_Writer), intent(inout) :: this
    end subroutine Null_Changes_Success_Writer_destroy

    subroutine Null_Changes_Success_Writer_write(this, i_step, changes_success)
        class(Null_Changes_Success_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Concrete_Changes_Success), intent(in) :: changes_success
    end subroutine Null_Changes_Success_Writer_write

!end implementation Null_Changes_Success_Writer

end module class_changes_writer
