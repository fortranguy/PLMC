module procedures_metropolis_algorithms_factory

use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use procedures_box_volume_change_factory, only: box_volume_change_create => create, &
    box_volume_change_destroy => destroy
use procedures_one_particle_move_factory, only: one_particle_translation_create => &
    create_translation, one_particle_rotation_create => create_rotation, &
    one_particle_move_destroy => destroy
use procedures_two_particles_switch_factory, only: two_particles_switch_create => create, &
    two_particles_switch_destroy => destroy
use procedures_one_particle_exchange_factory, only: one_particle_add_create => create_add, &
    one_particle_remove_create => create_remove, one_particle_exchange_destroy => destroy
use procedures_two_particles_transmutation_factory, only: two_particles_transmutation_create => &
    create, two_particles_transmutation_destroy => destroy
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithms_Wrapper

implicit none

private
public :: metropolis_algorithms_create, metropolis_algorithms_destroy, metropolis_algorithms_set

contains

    subroutine metropolis_algorithms_create(metropolis_algorithms, physical_model, changes)
        type(Metropolis_Algorithms_Wrapper), intent(out) :: metropolis_algorithms
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        call box_volume_change_create(metropolis_algorithms%box_volume_change, physical_model, &
            changes%changed_box_size)
        call one_particle_translation_create(metropolis_algorithms%one_particle_translation, &
            physical_model, changes%components)
        call one_particle_rotation_create(metropolis_algorithms%one_particle_rotation, &
            physical_model, changes%components)
        call two_particles_switch_create(metropolis_algorithms%two_particles_switch, physical_model)
        call one_particle_add_create(metropolis_algorithms%one_particle_add, physical_model, &
            changes)
        call one_particle_remove_create(metropolis_algorithms%one_particle_remove, physical_model, &
            changes)
        call two_particles_transmutation_create(metropolis_algorithms%two_particles_transmutation, &
            physical_model, changes)
    end subroutine metropolis_algorithms_create

    subroutine metropolis_algorithms_destroy(metropolis_algorithms)
        type(Metropolis_Algorithms_Wrapper), intent(inout) :: metropolis_algorithms

        call two_particles_transmutation_destroy(metropolis_algorithms%two_particles_transmutation)
        call one_particle_exchange_destroy(metropolis_algorithms%one_particle_remove)
        call one_particle_exchange_destroy(metropolis_algorithms%one_particle_add)
        call two_particles_switch_destroy(metropolis_algorithms%two_particles_switch)
        call one_particle_move_destroy(metropolis_algorithms%one_particle_rotation)
        call one_particle_move_destroy(metropolis_algorithms%one_particle_translation)
        call box_volume_change_destroy(metropolis_algorithms%box_volume_change)
    end subroutine metropolis_algorithms_destroy

    subroutine metropolis_algorithms_set(metropolis_algorithms)
        type(Metropolis_Algorithms_Wrapper), intent(inout) :: metropolis_algorithms

        call metropolis_algorithms%one_particle_translation%set_selector()
        call metropolis_algorithms%one_particle_rotation%set_selector()
        call metropolis_algorithms%two_particles_switch%set_selector()
        call metropolis_algorithms%one_particle_add%set_selector()
        call metropolis_algorithms%one_particle_remove%set_selector()
        call metropolis_algorithms%two_particles_transmutation%set_selector()
    end subroutine metropolis_algorithms_set

end module procedures_metropolis_algorithms_factory
