module procedures_metropolis_algorithms_factory

use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use procedures_one_particle_change_factory, only: one_particle_move_create => create_move, &
    one_particle_rotation_create => create_rotation, one_particle_change_destroy => destroy
use procedures_two_particles_switch_factory, only: two_particles_switch_create => create, &
    two_particles_switch_destroy => destroy, two_particles_switch_set => set
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithms_Wrapper

implicit none

private
public :: metropolis_algorithms_create, metropolis_algorithms_destroy, metropolis_algorithms_set

contains

    subroutine metropolis_algorithms_create(metropolis_algorithms, environment, mixture, &
        short_interactions, dipolar_interactions, change_components)
        type(Metropolis_Algorithms_Wrapper), intent(out) :: metropolis_algorithms
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions
        type(Changes_Component_Wrapper), intent(in) :: change_components(:)

        call one_particle_move_create(metropolis_algorithms%one_particle_move, &
            environment, mixture, short_interactions, dipolar_interactions, change_components)
        call one_particle_rotation_create(metropolis_algorithms%one_particle_rotation, &
            environment, mixture, short_interactions, dipolar_interactions, change_components)
        call two_particles_switch_create(metropolis_algorithms%two_particles_switch, environment, &
            mixture%components, short_interactions, dipolar_interactions)
    end subroutine metropolis_algorithms_create

    subroutine metropolis_algorithms_destroy(metropolis_algorithms)
        type(Metropolis_Algorithms_Wrapper), intent(inout) :: metropolis_algorithms

        call two_particles_switch_destroy(metropolis_algorithms%two_particles_switch)
        call one_particle_change_destroy(metropolis_algorithms%one_particle_rotation)
        call one_particle_change_destroy(metropolis_algorithms%one_particle_move)
    end subroutine metropolis_algorithms_destroy

    subroutine metropolis_algorithms_set(metropolis_algorithms, components)
        type(Metropolis_Algorithms_Wrapper), intent(inout) :: metropolis_algorithms
        type(Component_Wrapper), intent(in) :: components(:)

        call metropolis_algorithms%one_particle_move%set_selector()
        call metropolis_algorithms%one_particle_rotation%set_selector()

        call two_particles_switch_set(metropolis_algorithms%two_particles_switch, components)
    end subroutine metropolis_algorithms_set

end module procedures_metropolis_algorithms_factory
