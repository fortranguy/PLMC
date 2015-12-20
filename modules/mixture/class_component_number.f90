module class_component_number

use procedures_errors, only: error_exit

implicit none

private

    type, abstract, public :: Abstract_Component_Number
    private
        integer :: number
    contains
        procedure :: set => Abstract_set
        procedure :: get => Abstract_get
    end type Abstract_Component_Number

    type, extends(Abstract_Component_Number), public :: Concrete_Component_Number

    end type Concrete_Component_Number

    type, extends(Abstract_Component_Number), public :: Null_Component_Number
    contains
        procedure :: set => Null_set
        procedure :: get => Null_get
    end type Null_Component_Number

contains

!implementation Abstract_Component_Number

    subroutine Abstract_set(this, number)
        class(Abstract_Component_Number), intent(inout) :: this
        integer, intent(in) :: number

        if (number < 0) call error_exit("Abstract_Component_Number: number is negative.")
        this%number = number
    end subroutine Abstract_set

    pure function Abstract_get(this) result(number)
        class(Abstract_Component_Number), intent(in) :: this
        integer :: number

        number = this%number
    end function Abstract_get

!end implementation Abstract_Component_Number

!implementation Null_Component_Number

    subroutine Null_set(this, number)
        class(Null_Component_Number), intent(inout) :: this
        integer, intent(in) :: number
    end subroutine Null_set

    pure function Null_get(this) result(number)
        class(Null_Component_Number), intent(in) :: this
        integer :: number
        number = 0
    end function Null_get

!end implementation Null_Component_Number

end module class_component_number
