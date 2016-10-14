module procedures_exploring_writers_factory

use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_pair_potential_wrapper, only: Pair_Potential_Wrapper, Pair_Potentials_Line
use classes_volume_change_method, only: Abstract_Volume_Change_Method
use classes_particle_insertion_method, only: Abstract_Particle_Insertion_Method
use procedures_real_writer_factory, only: real_writer_create => create, &
    real_writer_destroy => destroy
use procedures_line_writer_factory, only: line_writer_create => create, &
    line_writer_destroy => destroy
use procedures_energies_writers_factory, only: energies_writers_create => create, &
    energies_writers_destroy => destroy
use types_exploring_writers_wrapper, only: Exploring_Writers_Wrapper
use procedures_exploration_inquirers, only: property_measure_pressure => measure_pressure
use procedures_exploration_inquirers, only: measure_chemical_potentials

implicit none

private
public :: create, destroy

contains

    !> @todo to update for multiple boxes
    subroutine create(writers, environment, wall_pairs, components, short_pairs, &
        volume_change_method, particle_insertion_method, visit_energies)
        type(Exploring_Writers_Wrapper), intent(out) :: writers
        type(Environment_Wrapper), intent(in) :: environment
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Pair_Potentials_Line), intent(in) :: short_pairs(:)
        class(Abstract_Volume_Change_Method), intent(in) :: volume_change_method
        class(Abstract_Particle_Insertion_Method), intent(in) :: particle_insertion_method
        logical, intent(in) :: visit_energies

        logical :: selector(size(components)), measure_pressure

        measure_pressure = property_measure_pressure(volume_change_method)
        call real_writer_create(writers%maximum_box_compression_delta, &
            "maximum_box_compression_delta.out", measure_pressure)
        call real_writer_create(writers%beta_pressure_excess, "beta_pressure_excess.out", &
            measure_pressure)
        call energies_writers_create(writers%energies, "", environment%external_field, wall_pairs, &
            components, short_pairs, visit_energies)
        selector = measure_chemical_potentials(particle_insertion_method)
        call line_writer_create(writers%insertion_successes, "insertion_successes.out", selector)
        call line_writer_create(writers%inv_pow_activities, "inv_pow_activities.out", selector)
    end subroutine create

    subroutine destroy(writers)
        type(Exploring_Writers_Wrapper), intent(inout) :: writers

        call energies_writers_destroy(writers%energies)
        call line_writer_destroy(writers%inv_pow_activities)
        call line_writer_destroy(writers%insertion_successes)
        call real_writer_destroy(writers%beta_pressure_excess)
        call real_writer_destroy(writers%maximum_box_compression_delta)
    end subroutine destroy

end module procedures_exploring_writers_factory
