module procedures_exploring_writers_factory

use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_total_moment_factory, only: set_are_dipolar
use types_pair_potential_wrapper, only: Pair_Potential_Wrapper, Pair_Potentials_Line
use classes_particle_insertion_method, only: Abstract_Particle_Insertion_Method
use procedures_real_writer_factory, only: real_writer_create => create, &
    real_writer_destroy => destroy
use procedures_line_writer_factory, only: line_writer_create => create, &
    line_writer_destroy => destroy
use procedures_triangle_writer_factory, only: triangle_writer_create => create, &
    triangle_writer_destroy => destroy
use types_exploring_writers_wrapper, only: Exploring_Writers_Wrapper
use procedures_property_inquirers, only: measure_chemical_potentials

implicit none

private
public :: create, destroy

contains

    subroutine create(writers, environment, wall_pairs, components, short_pairs, &
        particle_insertion_method)
        type(Exploring_Writers_Wrapper), intent(out) :: writers
        type(Environment_Wrapper), intent(in) :: environment
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Pair_Potentials_Line), intent(in) :: short_pairs(:)
        class(Abstract_Particle_Insertion_Method), intent(in) :: particle_insertion_method

        logical, dimension(size(components)) :: selector, are_dipolar

        selector = measure_chemical_potentials(particle_insertion_method)
        call line_writer_create(writers%insertion_successes, selector, "insertion_successes.out")
        call line_writer_create(writers%inv_pow_activities, selector, "inv_pow_activities.out")
        call line_writer_create(writers%field, environment%external_field, components, &
            "field_energies.out") !change filenames?
        call line_writer_create(writers%walls, wall_pairs, "walls_energies.out")
        call triangle_writer_create(writers%short_energies, short_pairs, "short_energies.out")
        call set_are_dipolar(are_dipolar, components)
        call triangle_writer_create(writers%dipolar_energies, are_dipolar, "dipolar_energies.out")
        call real_writer_create(writers%dipolar_mixture_energy, any(are_dipolar), &
            "dipolar_mixture_energy.out")
    end subroutine create

    subroutine destroy(writers)
        type(Exploring_Writers_Wrapper), intent(inout) :: writers

        call real_writer_destroy(writers%dipolar_mixture_energy)
        call triangle_writer_destroy(writers%dipolar_energies)
        call triangle_writer_destroy(writers%short_energies)
        call line_writer_destroy(writers%walls)
        call line_writer_destroy(writers%field)
        call line_writer_destroy(writers%inv_pow_activities)
        call line_writer_destroy(writers%insertion_successes)
    end subroutine destroy

end module procedures_exploring_writers_factory
