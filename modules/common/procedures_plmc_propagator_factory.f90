module procedures_plmc_propagator_factory

use types_generating_algorithms_wrapper, only: Generating_Algorithm_Pointer, &
    Generating_Algorithms_Wrapper
use classes_plmc_propagator, only: Abstract_PLMC_Propagator, Concrete_PLMC_Propagator, &
    Null_PLMC_Propagator

implicit none

private
public :: create, destroy

contains

    subroutine create(propagator, generating_algorithms)
        class(Abstract_PLMC_Propagator), allocatable, intent(out) :: propagator
        type(Generating_Algorithms_Wrapper), target, intent(in) :: generating_algorithms

        type(Generating_Algorithm_Pointer) :: algorithms(7)

        algorithms(1)%algorithm => generating_algorithms%box_volume_change
        algorithms(2)%algorithm => generating_algorithms%one_particle_translation
        algorithms(3)%algorithm => generating_algorithms%one_particle_rotation
        algorithms(4)%algorithm => generating_algorithms%two_particles_switch
        algorithms(5)%algorithm => generating_algorithms%one_particle_add
        algorithms(6)%algorithm => generating_algorithms%one_particle_remove
        algorithms(7)%algorithm => generating_algorithms%two_particles_transmutation
        allocate(Concrete_PLMC_Propagator :: propagator) !What about Null_?
        call propagator%construct(algorithms)
    end subroutine create

    subroutine destroy(propagator)
        class(Abstract_PLMC_Propagator), allocatable, intent(inout) :: propagator

        if (allocated(propagator)) then
            call propagator%destroy()
            deallocate(propagator)
        end if
    end subroutine destroy

end module procedures_plmc_propagator_factory
