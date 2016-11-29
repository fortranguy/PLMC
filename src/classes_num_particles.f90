module classes_num_particles

use procedures_errors, only: error_exit

implicit none

private

    type, abstract, public :: Abstract_Num_Particles
    private
        integer :: num_particles = 0
    contains
        procedure :: set => Abstract_set
        procedure :: get => Abstract_get
    end type Abstract_Num_Particles

    type, extends(Abstract_Num_Particles), public :: Concrete_Num_Particles

    end type Concrete_Num_Particles

    type, extends(Abstract_Num_Particles), public :: Null_Num_Particles
    contains
        procedure :: set => Null_set
        procedure :: get => Null_get
    end type Null_Num_Particles

contains

!implementation Abstract_Num_Particles

    subroutine Abstract_set(this, num_particles)
        class(Abstract_Num_Particles), intent(inout) :: this
        integer, intent(in) :: num_particles

        if (num_particles < 0) call error_exit("Abstract_Num_Particles: num_particles is negative.")
        this%num_particles = num_particles
    end subroutine Abstract_set

    pure function Abstract_get(this) result(num_particles)
        class(Abstract_Num_Particles), intent(in) :: this
        integer :: num_particles

        num_particles = this%num_particles
    end function Abstract_get

!end implementation Abstract_Num_Particles

!implementation Null_Num_Particles

    subroutine Null_set(this, num_particles)
        class(Null_Num_Particles), intent(inout) :: this
        integer, intent(in) :: num_particles
    end subroutine Null_set

    pure function Null_get(this) result(num_particles)
        class(Null_Num_Particles), intent(in) :: this
        integer :: num_particles
        num_particles = 0
    end function Null_get

!end implementation Null_Num_Particles

end module classes_num_particles
