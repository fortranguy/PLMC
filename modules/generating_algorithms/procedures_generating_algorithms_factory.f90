module procedures_generating_algorithms_factory

use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use procedures_box_volume_change_factory, only: box_volume_change_create => create, &
    box_volume_change_destroy => destroy
use procedures_one_particle_move_factory, only: one_particle_translation_create => &
    create_translation, one_particle_rotation_create => create_rotation, &
    one_particle_move_destroy => destroy
use procedures_one_particle_exchange_factory, only: one_particle_add_create => create_add, &
    one_particle_remove_create => create_remove, one_particle_exchange_destroy => destroy
use procedures_box_particles_swap_factory, only: box_particles_swap_identities_create => &
    create_identities, box_particles_swap_positions_create => create_positions, &
    box_particles_swap_destroy => destroy
use types_generating_algorithms_wrapper, only: Generating_Algorithms_Wrapper

implicit none

private
public :: generating_algorithms_create, generating_algorithms_destroy

contains

    subroutine generating_algorithms_create(generating_algorithms, physical_model, changes)
        type(Generating_Algorithms_Wrapper), intent(out) :: generating_algorithms
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        call box_volume_change_create(generating_algorithms%box_volume_change, physical_model, &
            changes%changed_box_size)
        call one_particle_translation_create(generating_algorithms%one_particle_translation, &
            physical_model, changes%gemc_components)
        call one_particle_rotation_create(generating_algorithms%one_particle_rotation, &
            physical_model, changes%gemc_components)
        call box_particles_swap_positions_create(generating_algorithms%two_particles_switch, &
            physical_model, changes)
        call one_particle_add_create(generating_algorithms%one_particle_add, physical_model, &
            changes)
        call one_particle_remove_create(generating_algorithms%one_particle_remove, physical_model, &
            changes)
        call box_particles_swap_identities_create(generating_algorithms%two_particles_transmutation, &
            physical_model, changes)
    end subroutine generating_algorithms_create

    subroutine generating_algorithms_destroy(generating_algorithms)
        type(Generating_Algorithms_Wrapper), intent(inout) :: generating_algorithms

        call box_particles_swap_destroy(generating_algorithms%two_particles_transmutation)
        call one_particle_exchange_destroy(generating_algorithms%one_particle_remove)
        call one_particle_exchange_destroy(generating_algorithms%one_particle_add)
        call box_particles_swap_destroy(generating_algorithms%two_particles_switch)
        call one_particle_move_destroy(generating_algorithms%one_particle_rotation)
        call one_particle_move_destroy(generating_algorithms%one_particle_translation)
        call box_volume_change_destroy(generating_algorithms%box_volume_change)
    end subroutine generating_algorithms_destroy

end module procedures_generating_algorithms_factory
