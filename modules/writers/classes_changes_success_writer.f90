module classes_changes_success_writer

use data_strings, only: max_line_length
use procedures_checks, only: check_string_not_empty
use classes_number_to_string, only: Abstract_Number_to_String, Concrete_Number_to_String, &
    Null_Number_to_String
use module_changes_success, only: Concrete_Changes_Success
use types_changes_success_writer_selector, only: Changes_Success_Writer_Selector

implicit none

private

    type, abstract, public :: Abstract_Changes_Success_Writer
    private
        class(Abstract_Number_to_String), allocatable :: string_translation
        class(Abstract_Number_to_String), allocatable :: string_rotation
        class(Abstract_Number_to_String), allocatable :: string_exchange
        integer :: file_unit = 0
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: write => Abstract_write
    end type Abstract_Changes_Success_Writer

    type, extends(Abstract_Changes_Success_Writer), public :: Concrete_Changes_Success_Writer

    end type Concrete_Changes_Success_Writer

    type, extends(Abstract_Changes_Success_Writer), public :: Null_Changes_Success_Writer
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: write => Null_write
    end type Null_Changes_Success_Writer

    type, public :: Changes_Success_Writer_Wrapper
        class(Abstract_Changes_Success_Writer), allocatable :: writer
    end type Changes_Success_Writer_Wrapper

contains

!implementation Abstract_Changes_Success_Writer

    subroutine Abstract_construct(this, filename, changes_selector)
        class(Abstract_Changes_Success_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        type(Changes_Success_Writer_Selector), intent(in) :: changes_selector

        character(len=:), allocatable :: legend

        call check_string_not_empty("Abstract_construct: filename", filename)
        open(newunit=this%file_unit, recl=max_line_length, file=filename, action="write")
        legend = "# i_step"
        if (changes_selector%write_translations) then
            allocate(Concrete_Number_to_String :: this%string_translation)
            legend = legend//"    translations"
        else
            allocate(Null_Number_to_String :: this%string_translation)
        end if
        if (changes_selector%write_rotations) then
            allocate(Concrete_Number_to_String :: this%string_rotation)
            legend = legend//"    rotations"
        else
            allocate(Null_Number_to_String :: this%string_rotation)
        end if
        if (changes_selector%write_exchanges) then
            allocate(Concrete_Number_to_String :: this%string_exchange)
            legend = legend//"    additions    removals"
        else
            allocate(Null_Number_to_String :: this%string_exchange)
        end if
        write(this%file_unit, *) legend
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Changes_Success_Writer), intent(inout) :: this

        if (allocated(this%string_exchange)) deallocate(this%string_exchange)
        if (allocated(this%string_rotation)) deallocate(this%string_rotation)
        if (allocated(this%string_translation)) deallocate(this%string_translation)
        close(this%file_unit)
    end subroutine Abstract_destroy

    subroutine Abstract_write(this, i_step, changes_success)
        class(Abstract_Changes_Success_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Concrete_Changes_Success), intent(in) :: changes_success

        write(this%file_unit, *) i_step, this%string_translation%get(changes_success%translation)//&
            this%string_rotation%get(changes_success%rotation)//&
            this%string_exchange%get(changes_success%add)//&
            this%string_exchange%get(changes_success%remove)
    end subroutine Abstract_write

!end implementation Abstract_Changes_Success_Writer

!implementation Null_Changes_Success_Writer

    subroutine Null_construct(this, filename, changes_selector)
        class(Null_Changes_Success_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        type(Changes_Success_Writer_Selector), intent(in) :: changes_selector
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Changes_Success_Writer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_write(this, i_step, changes_success)
        class(Null_Changes_Success_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Concrete_Changes_Success), intent(in) :: changes_success
    end subroutine Null_write

!end implementation Null_Changes_Success_Writer

end module classes_changes_success_writer
