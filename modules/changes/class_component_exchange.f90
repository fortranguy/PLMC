module class_component_exchange

use types_temporary_particle, only: Concrete_Temporary_Particle
use types_component_wrapper, only: Component_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Component_Exchange
    private
        type(Component_Wrapper), pointer :: component => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: add => Abstract_add
        procedure :: remove => Abstract_remove
    end type Abstract_Component_Exchange

    type, extends(Abstract_Component_Exchange), public :: Concrete_Component_Exchange

    end type Concrete_Component_Exchange

    type, extends(Abstract_Component_Exchange), public :: Null_Component_Exchange
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: add => Null_add
        procedure :: remove => Null_remove
    end type Null_Component_Exchange

contains

!implementation Abstract_Component_Exchange

    subroutine Abstract_construct(this, component)
        class(Abstract_Component_Exchange), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: component

        this%component => component
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Component_Exchange), intent(inout) :: this

        this%component => null()
    end subroutine Abstract_destroy

    subroutine Abstract_add(this, particle)
        class(Abstract_Component_Exchange), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle

        call this%component%number%set(this%component%number%get() + 1)
        call this%component%positions%add(particle%position)
        call this%component%orientations%add(particle%orientation)
    end subroutine Abstract_add

    subroutine Abstract_remove(this, i_particle)
        class(Abstract_Component_Exchange), intent(inout) :: this
        integer, intent(in) :: i_particle

        call this%component%orientations%remove(i_particle)
        call this%component%positions%remove(i_particle)
        call this%component%number%set(this%component%number%get() - 1)
    end subroutine Abstract_remove

!end implementation Abstract_Component_Exchange

!implementation Null_Component_Exchange

    subroutine Null_construct(this, component)
        class(Null_Component_Exchange), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: component
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Component_Exchange), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_add(this, particle)
        class(Null_Component_Exchange), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle
    end subroutine Null_add

    subroutine Null_remove(this, i_particle)
        class(Null_Component_Exchange), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_remove

!end implementation Null_Component_Exchange

end module class_component_exchange
