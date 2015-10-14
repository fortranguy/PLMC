module class_component_exchange

use types_temporary_particle, only: Concrete_Temporary_Particle
use types_component_wrapper, only: Component_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Component_Exchange
    private
        type(Component_Wrapper), pointer :: component => null()
    contains
        procedure :: construct => Abstract_Component_Exchange_construct
        procedure :: destroy => Abstract_Component_Exchange_destroy
        procedure :: add => Abstract_Component_Exchange_add
        procedure :: remove => Abstract_Component_Exchange_remove
    end type Abstract_Component_Exchange

    type, extends(Abstract_Component_Exchange), public :: Concrete_Component_Exchange

    end type Concrete_Component_Exchange

    type, extends(Abstract_Component_Exchange), public :: Null_Component_Exchange
    contains
        procedure :: construct => Null_Component_Exchange_construct
        procedure :: destroy => Null_Component_Exchange_destroy
        procedure :: add => Null_Component_Exchange_add
        procedure :: remove => Null_Component_Exchange_remove
    end type Null_Component_Exchange

contains

!implemetation Abstract_Component_Exchange

    subroutine Abstract_Component_Exchange_construct(this, component)
        class(Abstract_Component_Exchange), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: component

        this%component => component
    end subroutine Abstract_Component_Exchange_construct

    subroutine Abstract_Component_Exchange_destroy(this)
        class(Abstract_Component_Exchange), intent(inout) :: this

        this%component => null()
    end subroutine Abstract_Component_Exchange_destroy

    subroutine Abstract_Component_Exchange_add(this, particle)
        class(Abstract_Component_Exchange), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle

        call this%component%number%set(this%component%number%get() + 1)
        call this%component%positions%add(particle%position)
        call this%component%orientations%add(particle%orientation)
        call this%component%total_moment%add(this%component%dipolar_moments%get(&
            this%component%number%get()))
    end subroutine Abstract_Component_Exchange_add

    subroutine Abstract_Component_Exchange_remove(this, i_particle)
        class(Abstract_Component_Exchange), intent(inout) :: this
        integer, intent(in) :: i_particle

        call this%component%total_moment%remove(this%component%dipolar_moments%get(i_particle))
        call this%component%orientations%remove(i_particle)
        call this%component%positions%remove(i_particle)
        call this%component%number%set(this%component%number%get() - 1)
    end subroutine Abstract_Component_Exchange_remove

!end implemetation Abstract_Component_Exchange

!implemetation Null_Component_Exchange

    subroutine Null_Component_Exchange_construct(this, component)
        class(Null_Component_Exchange), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: component
    end subroutine Null_Component_Exchange_construct

    subroutine Null_Component_Exchange_destroy(this)
        class(Null_Component_Exchange), intent(inout) :: this
    end subroutine Null_Component_Exchange_destroy

    subroutine Null_Component_Exchange_add(this, particle)
        class(Null_Component_Exchange), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle
    end subroutine Null_Component_Exchange_add

    subroutine Null_Component_Exchange_remove(this, i_particle)
        class(Null_Component_Exchange), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_Component_Exchange_remove

!end implemetation Null_Component_Exchange

end module class_component_exchange
