module procedures_plmc_propagator_factory

use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only: tower_sampler_create => create, &
    tower_sampler_destroy => destroy
use types_component_wrapper, only: Component_Wrapper
use classes_generating_algorithm, only: Generating_Algorithm_Wrapper
use classes_plmc_propagator, only: Abstract_PLMC_Propagator, Concrete_PLMC_Propagator, &
    Null_PLMC_Propagator

implicit none

private
public :: create, destroy

contains

    !> @todo What about Null_?
    subroutine create(propagator, components, generating_algorithms)
        class(Abstract_PLMC_Propagator), allocatable, intent(out) :: propagator
        type(Component_Wrapper), intent(in) :: components(:, :)
        type(Generating_Algorithm_Wrapper), intent(in) :: generating_algorithms(:)

        class(Abstract_Tower_Sampler), allocatable :: selector

        allocate(Concrete_PLMC_Propagator :: propagator)
        call tower_sampler_create(selector, size(generating_algorithms), needed=.true.)
        call propagator%construct(components, generating_algorithms, selector)
        call tower_sampler_destroy(selector)
    end subroutine create

    subroutine destroy(propagator)
        class(Abstract_PLMC_Propagator), allocatable, intent(inout) :: propagator

        if (allocated(propagator)) then
            call propagator%destroy()
            deallocate(propagator)
        end if
    end subroutine destroy

end module procedures_plmc_propagator_factory
