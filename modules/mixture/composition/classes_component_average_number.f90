module classes_component_average_number

use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_num_particles, only: Abstract_Num_Particles
use classes_component_chemical_potential, only: Abstract_Component_Chemical_Potential

implicit none

private

    type, abstract, public :: Abstract_Component_Average_Number
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_get), deferred :: get
    end type Abstract_Component_Average_Number

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Component_Average_Number
            class(Abstract_Component_Average_Number), intent(inout) :: this
        end subroutine Abstract_destroy

        pure integer function Abstract_get(this) result(average_number)
        import :: Abstract_Component_Average_Number
            class(Abstract_Component_Average_Number), intent(in) :: this
        end function Abstract_get

    end interface

    type, extends(Abstract_Component_Average_Number), public :: Constant_Component_Average_Number
    private
        class(Abstract_Num_Particles), pointer :: num_particles => null()
    contains
        procedure :: construct => Constant_construct
        procedure :: destroy => Constant_destroy
        procedure :: get => Constant_get
    end type Constant_Component_Average_Number

    type, extends(Abstract_Component_Average_Number), public :: Variable_Component_Average_Number
    private
        class(Abstract_Parallelepiped_Domain), pointer :: accessible_domain => null()
        class(Abstract_Component_Chemical_Potential), pointer :: chemical_potential => null()
        integer :: average_num_particles
        integer :: accumulated_num_particles
        integer :: accumulation_period
    contains
        procedure :: construct => Variable_construct
        procedure :: destroy => Variable_destroy
        procedure :: get => Variable_get
    end type Variable_Component_Average_Number

    type, extends(Abstract_Component_Average_Number), public :: Null_Component_Average_Number
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: get => Null_get
    end type Null_Component_Average_Number

contains

!implementation Constant_Component_Average_Number

    subroutine Constant_construct(this, num_particles)
        class(Constant_Component_Average_Number), intent(out) :: this
        class(Abstract_Num_Particles), target, intent(in) :: num_particles

        this%num_particles => num_particles
    end subroutine Constant_construct

    subroutine Constant_destroy(this)
        class(Constant_Component_Average_Number), intent(inout) :: this

        this%num_particles => null()
    end subroutine Constant_destroy

    pure integer function Constant_get(this) result(average_number)
        class(Constant_Component_Average_Number), intent(in) :: this

        average_number = this%num_particles%get()
    end function Constant_get

!end implementation Constant_Component_Average_Number

!implementation Variable_Component_Average_Number

    subroutine Variable_construct(this, accessible_domain, chemical_potential)
        class(Variable_Component_Average_Number), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: accessible_domain
        class(Abstract_Component_Chemical_Potential), target, intent(in) :: chemical_potential

        this%accessible_domain => accessible_domain
        this%chemical_potential => chemical_potential
        this%accumulated_num_particles = 0
    end subroutine Variable_construct

    subroutine Variable_destroy(this)
        class(Variable_Component_Average_Number), intent(inout) :: this

        this%chemical_potential => null()
        this%accessible_domain => null()
    end subroutine Variable_destroy

    pure subroutine Variable_accumulate(this, i_step, num_particles)
        class(Variable_Component_Average_Number), intent(inout) :: this
        integer, intent(in) :: i_step
        integer, intent(in) :: num_particles

        this%accumulated_num_particles = this%accumulated_num_particles + num_particles
        if (mod(i_step, this%accumulation_period) == 0) then
            this%average_num_particles = (this%accumulated_num_particles + 1) / this%&
                accumulation_period
            this%accumulated_num_particles = 0
        end if
    end subroutine Variable_accumulate

    !> @warning: this value may fluctuate if the volume changes. Average volume needed?
    pure integer function Variable_get(this) result(average_number)
        class(Variable_Component_Average_Number), intent(in) :: this

        average_number = this%chemical_potential%get_density() * product(this%accessible_domain%&
            get_size())
    end function Variable_get

!end implementation Variable_Component_Average_Number

!implementation Null_Component_Average_Number

    subroutine Null_construct(this)
        class(Null_Component_Average_Number), intent(out) :: this
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Component_Average_Number), intent(inout) :: this
    end subroutine Null_destroy

    pure integer function Null_get(this) result(average_number)
        class(Null_Component_Average_Number), intent(in) :: this
        average_number = 0
    end function Null_get

!end implementation Null_Component_Average_Number

end module classes_component_average_number
