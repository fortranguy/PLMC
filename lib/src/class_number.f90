module class_number

use procedures_errors, only: error_exit

implicit none

private

    type, abstract, public :: Abstract_Number
    private
        integer :: number
    contains
        procedure :: set => Abstract_Number_set
        procedure :: get => Abstract_Number_get
    end type Abstract_Number

    type, extends(Abstract_Number), public :: Concrete_Number
        
    end type Concrete_Number

contains

    subroutine Abstract_Number_set(this, number)
        class(Abstract_Number), intent(inout) :: this
        integer, intent(in) :: number

        if (number < 0) call error_exit("Abstract_Number: number is negative.")
        this%number = number
    end subroutine Abstract_Number_set

    pure function Abstract_Number_get(this) result(number)
        class(Abstract_Number), intent(in) :: this
        integer :: number

        number = this%number
    end function Abstract_Number_get

end module class_number
