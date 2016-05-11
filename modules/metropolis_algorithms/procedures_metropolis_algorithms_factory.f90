module procedures_metropolis_algorithms_factory

use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_wrapper, only: Changes_Component_Wrapper
use procedures_one_particle_move_factory, only: one_particle_translation_create => &
    create_translation, one_particle_rotation_create => create_rotation, &
    one_particle_move_destroy => destroy
use procedures_two_particles_switch_factory, only: two_particles_switch_create => create, &
    two_particles_switch_destroy => destroy
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithms_Wrapper

implicit none

private
public :: metropolis_algorithms_create, metropolis_algorithms_destroy, metropolis_algorithms_set

contains

    subroutine metropolis_algorithms_create(metropolis_algorithms, physical_model, &
        change_components)
        type(Metropolis_Algorithms_Wrapper), intent(out) :: metropolis_algorithms
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Component_Wrapper), intent(in) :: change_components(:)

        call one_particle_translation_create(metropolis_algorithms%one_particle_translation, &
            physical_model, change_components)
        call one_particle_rotation_create(metropolis_algorithms%one_particle_rotation, &
            physical_model, change_components)
        call two_particles_switch_create(metropolis_algorithms%two_particles_switch, physical_model)
    end subroutine metropolis_algorithms_create

    subroutine metropolis_algorithms_destroy(metropolis_algorithms)
        type(Metropolis_Algorithms_Wrapper), intent(inout) :: metropolis_algorithms

        call two_particles_switch_destroy(metropolis_algorithms%two_particles_switch)
        call one_particle_move_destroy(metropolis_algorithms%one_particle_rotation)
        call one_particle_move_destroy(metropolis_algorithms%one_particle_translation)
    end subroutine metropolis_algorithms_destroy

    subroutine metropolis_algorithms_set(metropolis_algorithms)
        type(Metropolis_Algorithms_Wrapper), intent(inout) :: metropolis_algorithms

        call metropolis_algorithms%one_particle_translation%set_selector()
        call metropolis_algorithms%one_particle_rotation%set_selector()
        call metropolis_algorithms%two_particles_switch%set_selector()
    end subroutine metropolis_algorithms_set

end module procedures_metropolis_algorithms_factory
