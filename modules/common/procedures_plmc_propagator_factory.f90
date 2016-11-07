module procedures_plmc_propagator_factory

use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only: tower_sampler_create => create, &
    tower_sampler_destroy => destroy
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_factory, only: set_exist
use classes_generating_algorithm, only: Generating_Algorithm_Wrapper
use classes_plmc_propagator, only: Abstract_PLMC_Propagator, Concrete_PLMC_Propagator, &
    Null_PLMC_Propagator

implicit none

private
public :: create, destroy

contains

    !> @todo Is all(components_exist) too strict?
    subroutine create(propagator, components, generating_algorithms)
        class(Abstract_PLMC_Propagator), allocatable, intent(out) :: propagator
        type(Component_Wrapper), intent(in) :: components(:, :)
        type(Generating_Algorithm_Wrapper), intent(in) :: generating_algorithms(:)

        class(Abstract_Tower_Sampler), allocatable :: selector
        logical :: components_exist(size(components, 1), size(components, 2))

        call set_exist(components_exist, components)

        if (all(components_exist)) then
            allocate(Concrete_PLMC_Propagator :: propagator)
        else
            allocate(Null_PLMC_Propagator :: propagator)
        end if

        call tower_sampler_create(selector, size(generating_algorithms), all(components_exist))
        call propagator%construct(generating_algorithms, selector)
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
