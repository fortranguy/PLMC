module class_particles_exchange

use module_particles, only: Generic_Particle, Generic_Particles

implicit none

private

    type, public :: Particles_Exchange_Facade
    private
        type(Generic_Particles), pointer :: particles
    contains
        procedure :: construct => Particles_Exchange_Facade_construct
        procedure :: destroy => Particles_Exchange_Facade_destroy
        procedure :: add => Particles_Exchange_Facade_add
        procedure :: remove => Particles_Exchange_Facade_remove
    end type Particles_Exchange_Facade

contains

    subroutine Particles_Exchange_Facade_construct(this, particles)
        class(Particles_Exchange_Facade), intent(out) :: this
        type(Generic_Particles), target, intent(in) :: particles

        this%particles => particles
    end subroutine Particles_Exchange_Facade_construct

    subroutine Particles_Exchange_Facade_destroy(this)
        class(Particles_Exchange_Facade), intent(inout) :: this

        this%particles => null()
    end subroutine Particles_Exchange_Facade_destroy

    subroutine Particles_Exchange_Facade_add(this, particle)
        class(Particles_Exchange_Facade), intent(inout) :: this
        type(Generic_Particle), intent(in) :: particle

        call this%particles%number%set(this%particles%number%get() + 1)
        call this%particles%positions%add(particle%position)
        call this%particles%orientations%add(particle%orientation)
    end subroutine Particles_Exchange_Facade_add

    subroutine Particles_Exchange_Facade_remove(this, i_particle)
        class(Particles_Exchange_Facade), intent(inout) :: this
        integer, intent(in) :: i_particle

        call this%particles%orientations%remove(i_particle)
        call this%particles%positions%remove(i_particle)
        call this%particles%number%set(this%particles%number%get() - 1)
    end subroutine Particles_Exchange_Facade_remove

end module class_particles_exchange
