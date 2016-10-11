module procedures_one_particle_exchange_factory

use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only:tower_sampler_create => create, tower_sampler_destroy =>&
    destroy
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use classes_one_particle_exchange, only: Abstract_One_Particle_Exchange, Concrete_One_Particle_Add,&
    Concrete_One_Particle_Remove, Null_One_Particle_Exchange
use procedures_changes_factory, only: set_can_exchange

implicit none

private
public :: create_add, create_remove, destroy

contains

    subroutine create_add(one_particle_exchange, physical_model, changes)
        class(Abstract_One_Particle_Exchange), allocatable, intent(out) :: one_particle_exchange
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        class(Abstract_Tower_Sampler), allocatable :: selector
        logical :: can_exchange(size(physical_model%mixture%components))

        call set_can_exchange(can_exchange, physical_model%mixture%components)
        if (any(can_exchange)) then
            allocate(Concrete_One_Particle_Add :: one_particle_exchange)
        else
            allocate(Null_One_Particle_Exchange :: one_particle_exchange)
        end if

        call tower_sampler_create(selector, size(can_exchange), any(can_exchange))
        call one_particle_exchange%construct(physical_model%environment, physical_model%mixture, &
            physical_model%short_interactions, physical_model%dipolar_interactions_dynamic, &
            physical_model%dipolar_interactions_static, changes, can_exchange, selector)
        call tower_sampler_destroy(selector)
    end subroutine create_add

    subroutine create_remove(one_particle_exchange, physical_model, changes)
        class(Abstract_One_Particle_Exchange), allocatable, intent(out) :: one_particle_exchange
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        class(Abstract_Tower_Sampler), allocatable :: selector
        logical :: can_exchange(size(physical_model%mixture%components))

        call set_can_exchange(can_exchange, physical_model%mixture%components)
        if (any(can_exchange)) then
            allocate(Concrete_One_Particle_Remove :: one_particle_exchange)
        else
            allocate(Null_One_Particle_Exchange :: one_particle_exchange)
        end if

        call tower_sampler_create(selector, size(can_exchange), any(can_exchange))
        call one_particle_exchange%construct(physical_model%environment, physical_model%mixture, &
            physical_model%short_interactions, physical_model%dipolar_interactions_dynamic, &
            physical_model%dipolar_interactions_static, changes, can_exchange, selector)
        call tower_sampler_destroy(selector)
    end subroutine create_remove

    subroutine destroy(one_particle_exchange)
        class(Abstract_One_Particle_Exchange), allocatable, intent(inout) :: one_particle_exchange

        if (allocated(one_particle_exchange)) then
            call one_particle_exchange%destroy()
            deallocate(one_particle_exchange)
        end if
    end subroutine destroy

end module procedures_one_particle_exchange_factory
