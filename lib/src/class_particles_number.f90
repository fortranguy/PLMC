module class_particles_number

use module_error, only: error_exit

implicit none

private

    type, abstract, public :: Abstract_Particles_Number
    private
        integer :: num
    contains
        procedure :: set_num => Abstract_Particles_Number_set_num
        procedure :: get_num => Abstract_Particles_Number_get_num
    end type Abstract_Particles_Number

    type, extends(Abstract_Particles_Number), public :: Concrete_Particles_Number
        
    end type Concrete_Particles_Number

contains

    subroutine Abstract_Particles_Number_set_num(this, num)
        class(Abstract_Particles_Number), intent(out) :: this
        integer, intent(in) :: num

        if (num < 0) call error_exit("Particles Number is negative.")
        this%num = num
    end subroutine Abstract_Particles_Number_set_num

    pure function Abstract_Particles_Number_get_num(this) result(get_num)
        class(Abstract_Particles_Number), intent(in) :: this
        integer :: get_num

        get_num = this%num
    end function Abstract_Particles_Number_get_num

end module class_particles_number