module procedures_plmc_propagator_factory

use data_input_prefixes, only: changes_prefix
use json_module, only: json_file
use procedures_checks, only: check_data_found
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
    subroutine create(propagator, components, generating_algorithms, generating_data)
        class(Abstract_PLMC_Propagator), allocatable, intent(out) :: propagator
        type(Component_Wrapper), intent(in) :: components(:, :)
        type(Generating_Algorithm_Wrapper), intent(in) :: generating_algorithms(:)
        type(json_file), intent(inout) :: generating_data

        class(Abstract_Tower_Sampler), allocatable :: selector
        logical :: components_exist(size(components, 1), size(components, 2))
        integer :: tuning_period
        character(len=:), allocatable :: data_field
        logical :: data_found

        call set_exist(components_exist, components)

        if (all(components_exist)) then
            allocate(Concrete_PLMC_Propagator :: propagator)
            data_field = changes_prefix//"Propagator.tuning period"
            call generating_data%get(data_field, tuning_period, data_found)
            call check_data_found(data_field, data_found)
        else
            allocate(Null_PLMC_Propagator :: propagator)
            tuning_period = 0
        end if

        call tower_sampler_create(selector, size(generating_algorithms), all(components_exist))
        call propagator%construct(generating_algorithms, selector, tuning_period)
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
