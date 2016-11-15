module classes_average_num_particles

use procedures_checks, only: check_positive
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_num_particles, only: Abstract_Num_Particles
use classes_component_chemical_potential, only: Abstract_Component_Chemical_Potential

implicit none

private

    type, abstract, public :: Abstract_Average_Num_Particles
    private
        integer :: accumulation_period = 0
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_set), deferred :: set
        procedure(Abstract_accumulate), deferred :: accumulate
        procedure(Abstract_get), deferred :: get
        procedure, private :: set_accumulation_period => Abstract_set_accumulation_period
    end type Abstract_Average_Num_Particles

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Average_Num_Particles
            class(Abstract_Average_Num_Particles), intent(inout) :: this
        end subroutine Abstract_destroy

        pure subroutine Abstract_set(this)
        import :: Abstract_Average_Num_Particles
            class(Abstract_Average_Num_Particles), intent(inout) :: this
        end subroutine Abstract_set

        pure subroutine Abstract_accumulate(this, i_step)
        import :: Abstract_Average_Num_Particles
            class(Abstract_Average_Num_Particles), intent(inout) :: this
            integer, intent(in) :: i_step
        end subroutine Abstract_accumulate

        pure integer function Abstract_get(this) result(average_num_particles)
        import :: Abstract_Average_Num_Particles
            class(Abstract_Average_Num_Particles), intent(in) :: this
        end function Abstract_get

    end interface

    !> For canonical \( T, V, N \)  or isobaric \( T, p, N \) ensembles
    type, extends(Abstract_Average_Num_Particles), public :: Constant_Num_Particles
    private
        class(Abstract_Num_Particles), pointer :: num_particles => null()
    contains
        procedure :: construct => Constant_construct
        procedure :: destroy => Constant_destroy
        procedure :: set => Constant_set
        procedure :: accumulate => Constant_accumulate
        procedure :: get => Constant_get
    end type Constant_Num_Particles

    !> For grand canonical ensemble \( T, V, \mu \)
    type, extends(Abstract_Average_Num_Particles), public :: &
        Constant_Chemical_Potential_Num_Particles
    private
        class(Abstract_Parallelepiped_Domain), pointer :: accessible_domain => null()
        class(Abstract_Component_Chemical_Potential), pointer :: chemical_potential => null()
    contains
        procedure :: construct => Chemical_Potential_construct
        procedure :: destroy => Chemical_Potential_destroy
        procedure :: set => Chemical_Potential_set
        procedure :: accumulate => Chemical_Potential_accumulate
        procedure :: get => Chemical_Potential_get
    end type Constant_Chemical_Potential_Num_Particles

    type, extends(Abstract_Average_Num_Particles), public :: Concrete_Average_Num_Particles
    private
        class(Abstract_Num_Particles), pointer :: num_particles => null()
        integer :: average_num_particles = 0
        integer :: accumulated_num_particles = 0
    contains
        procedure :: construct => Concrete_construct
        procedure :: destroy => Concrete_destroy
        procedure :: set => Concrete_set
        procedure :: accumulate => Concrete_accumulate
        procedure :: get => Concrete_get
    end type Concrete_Average_Num_Particles

    type, extends(Abstract_Average_Num_Particles), public :: Null_Average_Num_Particles
    contains
        procedure :: destroy => Null_destroy
        procedure :: set => Null_set
        procedure :: accumulate => Null_accumulate
        procedure :: get => Null_get
    end type Null_Average_Num_Particles

contains

!implementation Abstract_Average_Num_Particles

    subroutine Abstract_set_accumulation_period(this, accumulation_period)
        class(Abstract_Average_Num_Particles), intent(inout) :: this
        integer, intent(in) :: accumulation_period

        call check_positive("Abstract_Average_Num_Particles", "accumulation_period", &
            accumulation_period)
        this%accumulation_period = accumulation_period
    end subroutine Abstract_set_accumulation_period

!end implementation Abstract_Average_Num_Particles

!implementation Constant_Num_Particles

    subroutine Constant_construct(this, num_particles, accumulation_period)
        class(Constant_Num_Particles), intent(out) :: this
        class(Abstract_Num_Particles), target, intent(in) :: num_particles
        integer, intent(in) :: accumulation_period

        this%num_particles => num_particles
        call this%set_accumulation_period(accumulation_period)
    end subroutine Constant_construct

    subroutine Constant_destroy(this)
        class(Constant_Num_Particles), intent(inout) :: this

        this%num_particles => null()
    end subroutine Constant_destroy

    pure subroutine Constant_set(this)
        class(Constant_Num_Particles), intent(inout) :: this

    end subroutine Constant_set

    pure subroutine Constant_accumulate(this, i_step)
        class(Constant_Num_Particles), intent(inout) :: this
        integer, intent(in) :: i_step

    end subroutine Constant_accumulate

    pure integer function Constant_get(this) result(average_num_particles)
        class(Constant_Num_Particles), intent(in) :: this

        average_num_particles = this%num_particles%get()
    end function Constant_get

!end implementation Constant_Num_Particles

!implementation Constant_Chemical_Potential_Num_Particles

    subroutine Chemical_Potential_construct(this, accessible_domain, chemical_potential, &
        accumulation_period)
        class(Constant_Chemical_Potential_Num_Particles), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: accessible_domain
        class(Abstract_Component_Chemical_Potential), target, intent(in) :: chemical_potential
        integer, intent(in) :: accumulation_period

        this%accessible_domain => accessible_domain
        this%chemical_potential => chemical_potential
        call this%set_accumulation_period(accumulation_period)
    end subroutine Chemical_Potential_construct

    subroutine Chemical_Potential_destroy(this)
        class(Constant_Chemical_Potential_Num_Particles), intent(inout) :: this

        this%chemical_potential => null()
        this%accessible_domain => null()
    end subroutine Chemical_Potential_destroy

    pure subroutine Chemical_Potential_set(this)
        class(Constant_Chemical_Potential_Num_Particles), intent(inout) :: this

    end subroutine Chemical_Potential_set

    pure subroutine Chemical_Potential_accumulate(this, i_step)
        class(Constant_Chemical_Potential_Num_Particles), intent(inout) :: this
        integer, intent(in) :: i_step

    end subroutine Chemical_Potential_accumulate

    !> @warning What if the field changes the number of particles?
    pure integer function Chemical_Potential_get(this) result(average_num_particles)
        class(Constant_Chemical_Potential_Num_Particles), intent(in) :: this

        average_num_particles = ceiling(this%chemical_potential%get_density() * &
            product(this%accessible_domain%get_size()))
    end function Chemical_Potential_get

!end implementation Constant_Chemical_Potential_Num_Particles

!implementation Concrete_Average_Num_Particles

    subroutine Concrete_construct(this, num_particles, accumulation_period)
        class(Concrete_Average_Num_Particles), intent(out) :: this
        class(Abstract_Num_Particles), target, intent(in) :: num_particles
        integer, intent(in) :: accumulation_period

        this%num_particles => num_particles
        call this%set_accumulation_period(accumulation_period)
    end subroutine Concrete_construct

    subroutine Concrete_destroy(this)
        class(Concrete_Average_Num_Particles), intent(inout) :: this

        this%num_particles => null()
    end subroutine Concrete_destroy

    pure subroutine Concrete_set(this)
        class(Concrete_Average_Num_Particles), intent(inout) :: this

        this%average_num_particles = this%num_particles%get()
    end subroutine Concrete_set

    !> @note minimum of this%average_num_particles is 1. Is it useless?
    pure subroutine Concrete_accumulate(this, i_step)
        class(Concrete_Average_Num_Particles), intent(inout) :: this
        integer, intent(in) :: i_step

        this%accumulated_num_particles = this%accumulated_num_particles + this%num_particles%get()
        if (mod(i_step, this%accumulation_period) == 0) then
            this%average_num_particles = this%accumulated_num_particles / this%accumulation_period
            if (this%average_num_particles < 1) this%average_num_particles = 1
            this%accumulated_num_particles = 0
        end if
    end subroutine Concrete_accumulate

    pure integer function Concrete_get(this) result(average_num_particles)
        class(Concrete_Average_Num_Particles), intent(in) :: this

        average_num_particles = this%average_num_particles
    end function Concrete_get

!end implementation Concrete_Average_Num_Particles

!implementation Null_Average_Num_Particles

    subroutine Null_destroy(this)
        class(Null_Average_Num_Particles), intent(inout) :: this
    end subroutine Null_destroy

    pure subroutine Null_set(this)
        class(Null_Average_Num_Particles), intent(inout) :: this
    end subroutine Null_set

    pure subroutine Null_accumulate(this, i_step)
        class(Null_Average_Num_Particles), intent(inout) :: this
        integer, intent(in) :: i_step
    end subroutine Null_accumulate

    pure integer function Null_get(this) result(average_num_particles)
        class(Null_Average_Num_Particles), intent(in) :: this
        average_num_particles = 0
    end function Null_get

!end implementation Null_Average_Num_Particles

end module classes_average_num_particles
