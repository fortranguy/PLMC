module procedures_two_particles_transmutation_factory

use classes_hetero_couples, only: Abstract_Hetero_Couples
use procedures_hetero_couples_factory, only: hetero_couples_create_full => create_full, &
    hetero_couples_destroy => destroy
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only:tower_sampler_create => create, tower_sampler_destroy =>&
    destroy
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use classes_two_particles_transmutation, only: Abstract_Two_Particles_Transmutation, &
    Concrete_Two_Particles_Transmutation, Null_Two_Particles_Transmutation
use procedures_changes_factory, only: set_can_exchange

implicit none

private
public :: create, destroy

contains

    subroutine create(two_particles_transmutation, physical_model, changes)
        class(Abstract_Two_Particles_Transmutation), allocatable, intent(out) :: &
            two_particles_transmutation
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        class(Abstract_Hetero_Couples), allocatable :: couples
        class(Abstract_Tower_Sampler), allocatable :: selector
        logical :: can_exchange(size(physical_model%mixture%gemc_components, 1))

        call set_can_exchange(can_exchange, physical_model%mixture%gemc_components(:, 1))
        if (count(can_exchange) > 1) then
            allocate(Concrete_Two_Particles_Transmutation :: two_particles_transmutation)
        else
            allocate(Null_Two_Particles_Transmutation :: two_particles_transmutation)
        end if

        call hetero_couples_create_full(couples, size(can_exchange))
        call tower_sampler_create(selector, couples%get_num(), count(can_exchange) > 1)
        call two_particles_transmutation%construct(physical_model%environment, physical_model%&
            mixture, physical_model%short_interactions, physical_model%&
            dipolar_interactions_dynamic, physical_model%dipolar_interactions_static, changes, &
            can_exchange, couples, selector)
        call tower_sampler_destroy(selector)
        call hetero_couples_destroy(couples)
    end subroutine create

    subroutine destroy(two_particles_transmutation)
        class(Abstract_Two_Particles_Transmutation), allocatable, intent(inout) :: &
            two_particles_transmutation

        if (allocated(two_particles_transmutation)) then
            call two_particles_transmutation%destroy()
            deallocate(two_particles_transmutation)
        end if
    end subroutine destroy

end module procedures_two_particles_transmutation_factory
