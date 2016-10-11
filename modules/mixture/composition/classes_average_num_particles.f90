module classes_average_num_particles

use procedures_checks, only: check_positive
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_num_particles, only: Abstract_Num_Particles
use classes_component_chemical_potential, only: Abstract_Component_Chemical_Potential

implicit none

private

    type, abstract, public :: Abstract_Average_Num_Particles
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_accumulate), deferred :: accumulate
        procedure(Abstract_get), deferred :: get
    end type Abstract_Average_Num_Particles

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Average_Num_Particles
            class(Abstract_Average_Num_Particles), intent(inout) :: this
        end subroutine Abstract_destroy

        pure subroutine Abstract_accumulate(this, i_step, num_particles)
        import :: Abstract_Average_Num_Particles
            class(Abstract_Average_Num_Particles), intent(inout) :: this
            integer, intent(in) :: i_step
            integer, intent(in) :: num_particles
        end subroutine Abstract_accumulate

        pure integer function Abstract_get(this) result(num_particles)
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
        procedure :: accumulate => Constant_accumulate
        procedure :: get => Constant_get
    end type Constant_Num_Particles

    !> For grand canonical ensemble \( T, V, \mu \)
    !> @todo
    !> What if the field is on and the number of particles changes?
    type, extends(Abstract_Average_Num_Particles), public :: &
        Constant_Chemical_Potential_Num_Particles
    private
        class(Abstract_Parallelepiped_Domain), pointer :: accessible_domain => null()
        class(Abstract_Component_Chemical_Potential), pointer :: chemical_potential => null()
    contains
        procedure :: construct => Chemical_Potential_construct
        procedure :: destroy => Chemical_Potential_destroy
        procedure :: accumulate => Chemical_Potential_accumulate
        procedure :: get => Chemical_Potential_get
    end type Constant_Chemical_Potential_Num_Particles

    !> @todo
    !> better initial this%num_particles?
    type, extends(Abstract_Average_Num_Particles), public :: Concrete_Average_Num_Particles
    private
        integer :: num_particles = 1
        integer :: accumulated_num_particles = 0
        integer :: accumulation_period = 0
    contains
        procedure :: construct => Concrete_construct
        procedure :: destroy => Concrete_destroy
        procedure :: accumulate => Concrete_accumulate
        procedure :: get => Concrete_get
    end type Concrete_Average_Num_Particles

    type, extends(Abstract_Average_Num_Particles), public :: Null_Average_Num_Particles
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: accumulate => Null_accumulate
        procedure :: get => Null_get
    end type Null_Average_Num_Particles

contains

!implementation Constant_Num_Particles

    subroutine Constant_construct(this, num_particles)
        class(Constant_Num_Particles), intent(out) :: this
        class(Abstract_Num_Particles), target, intent(in) :: num_particles

        this%num_particles => num_particles
    end subroutine Constant_construct

    subroutine Constant_destroy(this)
        class(Constant_Num_Particles), intent(inout) :: this

        this%num_particles => null()
    end subroutine Constant_destroy

    pure subroutine Constant_accumulate(this, i_step, num_particles)
        class(Constant_Num_Particles), intent(inout) :: this
        integer, intent(in) :: i_step
        integer, intent(in) :: num_particles

    end subroutine Constant_accumulate

    pure integer function Constant_get(this) result(num_particles)
        class(Constant_Num_Particles), intent(in) :: this

        num_particles = this%num_particles%get()
    end function Constant_get

!end implementation Constant_Num_Particles

!implementation Constant_Chemical_Potential_Num_Particles

    subroutine Chemical_Potential_construct(this, accessible_domain, chemical_potential)
        class(Constant_Chemical_Potential_Num_Particles), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: accessible_domain
        class(Abstract_Component_Chemical_Potential), target, intent(in) :: chemical_potential

        this%accessible_domain => accessible_domain
        this%chemical_potential => chemical_potential
    end subroutine Chemical_Potential_construct

    subroutine Chemical_Potential_destroy(this)
        class(Constant_Chemical_Potential_Num_Particles), intent(inout) :: this

        this%chemical_potential => null()
        this%accessible_domain => null()
    end subroutine Chemical_Potential_destroy

    pure subroutine Chemical_Potential_accumulate(this, i_step, num_particles)
        class(Constant_Chemical_Potential_Num_Particles), intent(inout) :: this
        integer, intent(in) :: i_step
        integer, intent(in) :: num_particles

    end subroutine Chemical_Potential_accumulate

    pure integer function Chemical_Potential_get(this) result(num_particles)
        class(Constant_Chemical_Potential_Num_Particles), intent(in) :: this

        num_particles = this%chemical_potential%get_density() * &
            product(this%accessible_domain%get_size())
    end function Chemical_Potential_get

!end implementation Constant_Chemical_Potential_Num_Particles

!implementation Concrete_Average_Num_Particles

    subroutine Concrete_construct(this, accumulation_period)
        class(Concrete_Average_Num_Particles), intent(out) :: this
        integer, intent(in) :: accumulation_period

        call check_positive("Concrete_Average_Num_Particles", "accumulation_period", &
            accumulation_period)
        this%accumulation_period = accumulation_period
    end subroutine Concrete_construct

    subroutine Concrete_destroy(this)
        class(Concrete_Average_Num_Particles), intent(inout) :: this

    end subroutine Concrete_destroy

    pure subroutine Concrete_accumulate(this, i_step, num_particles)
        class(Concrete_Average_Num_Particles), intent(inout) :: this
        integer, intent(in) :: i_step
        integer, intent(in) :: num_particles

        this%accumulated_num_particles = this%accumulated_num_particles + num_particles
        if (mod(i_step, this%accumulation_period) == 0) then
            this%num_particles = this%accumulated_num_particles / this%accumulation_period
            if (this%num_particles < 1) this%num_particles = 1
            this%accumulated_num_particles = 0
        end if
    end subroutine Concrete_accumulate

    pure integer function Concrete_get(this) result(num_particles)
        class(Concrete_Average_Num_Particles), intent(in) :: this

        num_particles = this%num_particles
    end function Concrete_get

!end implementation Concrete_Average_Num_Particles

!implementation Null_Average_Num_Particles

    subroutine Null_construct(this)
        class(Null_Average_Num_Particles), intent(out) :: this
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Average_Num_Particles), intent(inout) :: this
    end subroutine Null_destroy

    pure subroutine Null_accumulate(this, i_step, num_particles)
        class(Null_Average_Num_Particles), intent(inout) :: this
        integer, intent(in) :: i_step
        integer, intent(in) :: num_particles
    end subroutine Null_accumulate

    pure integer function Null_get(this) result(num_particles)
        class(Null_Average_Num_Particles), intent(in) :: this
        num_particles = 0
    end function Null_get

!end implementation Null_Average_Num_Particles

end module classes_average_num_particles
