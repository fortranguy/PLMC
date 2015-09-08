module class_particles_number

use procedures_errors, only: error_exit

implicit none

private

    type, abstract, public :: Abstract_Particles_Number
    private
        integer :: num
    contains
        procedure :: set => Abstract_Particles_Number_set
        procedure :: get => Abstract_Particles_Number_get
    end type Abstract_Particles_Number
    
    type, public :: Abstract_Particles_Number_Pointer
        class(Abstract_Particles_Number), pointer :: ptr
    end type Abstract_Particles_Number_Pointer

    type, extends(Abstract_Particles_Number), public :: Concrete_Particles_Number
        
    end type Concrete_Particles_Number

contains

    subroutine Abstract_Particles_Number_set(this, num)
        class(Abstract_Particles_Number), intent(out) :: this
        integer, intent(in) :: num

        if (num < 0) call error_exit("Abstract_Particles_Number: num is negative.")
        this%num = num
    end subroutine Abstract_Particles_Number_set

    pure function Abstract_Particles_Number_get(this) result(num)
        class(Abstract_Particles_Number), intent(in) :: this
        integer :: num

        num = this%num
    end function Abstract_Particles_Number_get

end module class_particles_number
