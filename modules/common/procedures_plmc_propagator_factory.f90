module procedures_plmc_propagator_factory

use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithm_Pointer, &
    Metropolis_Algorithms_Wrapper
use classes_plmc_propagator, only: Abstract_PLMC_Propagator, Concrete_PLMC_Propagator, &
    Null_PLMC_Propagator

implicit none

private
public :: create, destroy

contains

    subroutine create(propagator, metropolis_algorithms)
        class(Abstract_PLMC_Propagator), allocatable, intent(out) :: propagator
        type(Metropolis_Algorithms_Wrapper), target, intent(in) :: metropolis_algorithms

        type(Metropolis_Algorithm_Pointer) :: algorithms(3)

        algorithms(1)%algorithm => metropolis_algorithms%one_particle_translation
        algorithms(2)%algorithm => metropolis_algorithms%one_particle_rotation
        algorithms(3)%algorithm => metropolis_algorithms%two_particles_switch

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
