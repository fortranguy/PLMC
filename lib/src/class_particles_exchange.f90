module class_particles_exchange

use types_temporary_particle, only: Concrete_Temporary_Particle
use types_particles_wrapper, only: Particles_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Particles_Exchange
    private
        type(Particles_Wrapper), pointer :: particles => null()
    contains
        procedure :: construct => Abstract_Particles_Exchange_construct
        procedure :: destroy => Abstract_Particles_Exchange_destroy
        procedure :: add => Abstract_Particles_Exchange_add
        procedure :: remove => Abstract_Particles_Exchange_remove
    end type Abstract_Particles_Exchange

    type, extends(Abstract_Particles_Exchange), public :: Concrete_Particles_Exchange

    end type Concrete_Particles_Exchange

    type, extends(Abstract_Particles_Exchange), public :: Null_Particles_Exchange
    contains
        procedure :: construct => Null_Particles_Exchange_construct
        procedure :: destroy => Null_Particles_Exchange_destroy
        procedure :: add => Null_Particles_Exchange_add
        procedure :: remove => Null_Particles_Exchange_remove
    end type Null_Particles_Exchange

contains

!implemetation Abstract_Particles_Exchange

    subroutine Abstract_Particles_Exchange_construct(this, particles)
        class(Abstract_Particles_Exchange), intent(out) :: this
        type(Particles_Wrapper), target, intent(in) :: particles

        this%particles => particles
    end subroutine Abstract_Particles_Exchange_construct

    subroutine Abstract_Particles_Exchange_destroy(this)
        class(Abstract_Particles_Exchange), intent(inout) :: this

        this%particles => null()
    end subroutine Abstract_Particles_Exchange_destroy

    subroutine Abstract_Particles_Exchange_add(this, particle)
        class(Abstract_Particles_Exchange), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle

        call this%particles%number%set(this%particles%number%get() + 1)
        call this%particles%positions%add(particle%position)
        call this%particles%orientations%add(particle%orientation)
        call this%particles%total_moment%add(this%particles%dipolar_moments%get(&
            this%particles%number%get()))
    end subroutine Abstract_Particles_Exchange_add

    subroutine Abstract_Particles_Exchange_remove(this, i_particle)
        class(Abstract_Particles_Exchange), intent(inout) :: this
        integer, intent(in) :: i_particle

        call this%particles%total_moment%remove(this%particles%dipolar_moments%get(i_particle))
        call this%particles%orientations%remove(i_particle)
        call this%particles%positions%remove(i_particle)
        call this%particles%number%set(this%particles%number%get() - 1)
    end subroutine Abstract_Particles_Exchange_remove

!end implemetation Abstract_Particles_Exchange

!implemetation Null_Particles_Exchange

    subroutine Null_Particles_Exchange_construct(this, particles)
        class(Null_Particles_Exchange), intent(out) :: this
        type(Particles_Wrapper), target, intent(in) :: particles
    end subroutine Null_Particles_Exchange_construct

    subroutine Null_Particles_Exchange_destroy(this)
        class(Null_Particles_Exchange), intent(inout) :: this
    end subroutine Null_Particles_Exchange_destroy

    subroutine Null_Particles_Exchange_add(this, particle)
        class(Null_Particles_Exchange), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle
    end subroutine Null_Particles_Exchange_add

    subroutine Null_Particles_Exchange_remove(this, i_particle)
        class(Null_Particles_Exchange), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_Particles_Exchange_remove

!end implemetation Null_Particles_Exchange

end module class_particles_exchange
