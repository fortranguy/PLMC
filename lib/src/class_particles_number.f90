module class_particles_number

use procedures_errors, only: error_exit

implicit none

private

    type, abstract, public :: Abstract_Particles_Number
    private
        integer :: number
    contains
        procedure :: set => Abstract_Particles_Number_set
        procedure :: get => Abstract_Particles_Number_get
    end type Abstract_Particles_Number

    type, extends(Abstract_Particles_Number), public :: Concrete_Particles_Number

    end type Concrete_Particles_Number

    type, extends(Abstract_Particles_Number), public :: Null_Particles_Number
    contains
        procedure :: set => Null_Particles_Number_set
        procedure :: get => Null_Particles_Number_get
    end type Null_Particles_Number

contains

!implementation Abstract_Particles_Number

    subroutine Abstract_Particles_Number_set(this, number)
        class(Abstract_Particles_Number), intent(inout) :: this
        integer, intent(in) :: number

        if (number < 0) call error_exit("Abstract_Particles_Number: number is negative.")
        this%number = number
    end subroutine Abstract_Particles_Number_set

    pure function Abstract_Particles_Number_get(this) result(number)
        class(Abstract_Particles_Number), intent(in) :: this
        integer :: number

        number = this%number
    end function Abstract_Particles_Number_get

!end implementation Abstract_Particles_Number

!implementation Null_Particles_Number

    subroutine Null_Particles_Number_set(this, number)
        class(Null_Particles_Number), intent(inout) :: this
        integer, intent(in) :: number
    end subroutine Null_Particles_Number_set

    pure function Null_Particles_Number_get(this) result(number)
        class(Null_Particles_Number), intent(in) :: this
        integer :: number
        number = 0
    end function Null_Particles_Number_get

!end implementation Null_Particles_Number

end module class_particles_number
