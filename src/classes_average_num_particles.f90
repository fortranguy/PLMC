module classes_average_num_particles

use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_num_particles, only: Abstract_Num_Particles
use classes_component_chemical_potential, only: Abstract_Component_Chemical_Potential
use types_component_wrapper, only: Component_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Average_Num_Particles
    private
        integer :: average_num_particles = 0
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_set), deferred :: set
        procedure :: get => Abstract_get
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

    end interface

    !> For canonical \( T, V, N \)  or isobaric \( T, p, N \) ensembles
    type, extends(Abstract_Average_Num_Particles), public :: Constant_Num_Particles
    private
        class(Abstract_Num_Particles), pointer :: num_particles => null()
    contains
        procedure :: construct => Constant_construct
        procedure :: destroy => Constant_destroy
        procedure :: set => Constant_set
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
    end type Constant_Chemical_Potential_Num_Particles

    !> @note For GEMC, averaging the number of particles during thermalisation seems to affect
    !> the Markov chain too much. A temporary solution is to equally divide the number of particles
    !> among components.
    type, extends(Abstract_Average_Num_Particles), public :: Equipartition_Num_Particles
    private
        type(Component_Wrapper), pointer :: components(:) => null()
    contains
        procedure :: construct => Equipartition_construct
        procedure :: destroy => Equipartition_destroy
        procedure :: set => Equipartition_set
    end type Equipartition_Num_Particles

    type, extends(Abstract_Average_Num_Particles), public :: Null_Average_Num_Particles
    contains
        procedure :: destroy => Null_destroy
        procedure :: set => Null_set
    end type Null_Average_Num_Particles

contains

!implementation Abstract_Average_Num_Particles

   pure integer function Abstract_get(this) result(average_num_particles)
        class(Abstract_Average_Num_Particles), intent(in) :: this

        average_num_particles = this%average_num_particles
    end function Abstract_get

!end implementation Abstract_Average_Num_Particles

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

    pure subroutine Constant_set(this)
        class(Constant_Num_Particles), intent(inout) :: this

        this%average_num_particles = this%num_particles%get()
    end subroutine Constant_set

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

    !> @warning What if the field changes the number of particles
    !> or the accessible domain size changes?
    pure subroutine Chemical_Potential_set(this)
        class(Constant_Chemical_Potential_Num_Particles), intent(inout) :: this

        this%average_num_particles = ceiling(this%chemical_potential%get_density() * &
            product(this%accessible_domain%get_size()))
    end subroutine Chemical_Potential_set

!end implementation Constant_Chemical_Potential_Num_Particles

!implementation Equipartition_Num_Particles

    subroutine Equipartition_construct(this, componens)
        class(Equipartition_Num_Particles), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: componens(:)

        this%components => componens
    end subroutine Equipartition_construct

    subroutine Equipartition_destroy(this)
        class(Equipartition_Num_Particles), intent(inout) :: this

        this%components => null()
    end subroutine Equipartition_destroy

    !> @todo What if average_num_particles == 0?
    !> cf. [[procedures_selectors_resetters]]
    pure subroutine Equipartition_set(this)
        class(Equipartition_Num_Particles), intent(inout) :: this

        integer :: i_component, num_particles

        num_particles = 0
        do i_component = 1, size(this%components)
            num_particles = num_particles + this%components(i_component)%num_particles%get()
        end do
        this%average_num_particles = num_particles / size(this%components)
        if (this%average_num_particles == 0) this%average_num_particles = 1
    end subroutine Equipartition_set

!end implementation Equipartition_Num_Particles

!implementation Null_Average_Num_Particles

    subroutine Null_destroy(this)
        class(Null_Average_Num_Particles), intent(inout) :: this
    end subroutine Null_destroy

    pure subroutine Null_set(this)
        class(Null_Average_Num_Particles), intent(inout) :: this
    end subroutine Null_set

!end implementation Null_Average_Num_Particles

end module classes_average_num_particles
