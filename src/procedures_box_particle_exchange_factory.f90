module procedures_box_particle_exchange_factory

use procedures_errors, only: error_exit
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only:tower_sampler_create => create, tower_sampler_destroy =>&
    destroy
use procedures_mixture_properties, only: set_can_exchange
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use classes_generating_algorithm, only: Abstract_Generating_Algorithm, Null_Generating_Algorithm
use classes_box_particle_exchange, only: Box_Particle_Add, Box_Particle_Remove

implicit none

private
public :: create_add, create_remove

contains

    subroutine create_add(particle_add, physical_model, changes)
        class(Abstract_Generating_Algorithm), allocatable, intent(out) :: particle_add
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        class(Abstract_Tower_Sampler), allocatable :: selectors(:)
        logical :: can_exchange(size(physical_model%mixture%components, 1), &
            size(physical_model%mixture%components, 2))

        call set_can_exchange(can_exchange, physical_model%mixture%components)
        if (any(can_exchange)) then
            allocate(Box_Particle_Add :: particle_add)
        else
            allocate(Null_Generating_Algorithm :: particle_add)
        end if

        call tower_sampler_create(selectors, size(can_exchange, 2), size(can_exchange, 1), &
            any(can_exchange))
        select type (particle_add)
            type is (Box_Particle_Add)
                call particle_add%construct(physical_model%environment, physical_model%mixture, &
                    physical_model%short_interactions, physical_model%dipolar_interactions_dynamic,&
                    physical_model%dipolar_interactions_static, changes, can_exchange, &
                    selectors)
            type is (Null_Generating_Algorithm)
            class default
                call error_exit("procedures_box_particle_exchange_factory: create_add: "//&
                    "particle_add: unknown type.")
        end select
        call tower_sampler_destroy(selectors)
    end subroutine create_add

    subroutine create_remove(particle_remove, physical_model, changes)
        class(Abstract_Generating_Algorithm), allocatable, intent(out) :: particle_remove
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        class(Abstract_Tower_Sampler), allocatable :: selectors(:)
        logical :: can_exchange(size(physical_model%mixture%components, 1), &
            size(physical_model%mixture%components, 2))

        call set_can_exchange(can_exchange, physical_model%mixture%components)
        if (any(can_exchange)) then
            allocate(Box_Particle_Remove :: particle_remove)
        else
            allocate(Null_Generating_Algorithm :: particle_remove)
        end if

        call tower_sampler_create(selectors, size(can_exchange, 2), size(can_exchange, 1), &
            any(can_exchange))
        select type (particle_remove)
            type is (Box_Particle_Remove)
                call particle_remove%construct(physical_model%environment, physical_model%mixture, &
                    physical_model%short_interactions, physical_model%dipolar_interactions_dynamic,&
                    physical_model%dipolar_interactions_static, changes, can_exchange, &
                    selectors)
            type is (Null_Generating_Algorithm)
            class default
                call error_exit("procedures_box_particle_exchange_factory: create_remove: "//&
                    "particle_remove: unknown type.")
        end select
        call tower_sampler_destroy(selectors)
    end subroutine create_remove

end module procedures_box_particle_exchange_factory
