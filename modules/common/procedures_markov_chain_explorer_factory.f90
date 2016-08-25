module procedures_markov_chain_explorer_factory

use data_input_prefixes, only: particle_insertion_prefix, volume_change_prefix
use json_module, only: json_file
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use procedures_box_factory, only: box_create => create, box_destroy => destroy
use procedures_mixture_factory, only: set_have_positions, set_have_orientations
use procedures_random_coordinates_factory, only: random_coordinates_create => create, &
    random_coordinates_destroy => destroy
use classes_maximum_box_compression, only: Abstract_Maximum_Box_Compression
use procedures_maximum_box_compression_factory, only: maximum_box_compression_create => create, &
    maximum_box_compression_destroy => destroy
use procedures_maximum_box_compression_explorer_factory, only: &
    maximum_box_compression_explorer_create => create, &
    maximum_box_compression_explorer_destroy => destroy
use procedures_volume_change_method_factory, only: volume_change_method_create => create, &
    volume_change_method_destroy => destroy
use procedures_particle_insertion_method_factory, only: particle_insertion_method_create => create,&
    particle_insertion_method_destroy => destroy
use types_markov_chain_explorer_wrapper, only: Markov_Chain_Explorer_Wrapper
use procedures_property_inquirers, only: measure_maximum_compression, measure_pressure, &
    measure_chemical_potentials

implicit none

private
public :: create, destroy

contains

    subroutine create(markov_chain_explorer, physical_model, exploring_data)
        type(Markov_Chain_Explorer_Wrapper), intent(out) :: markov_chain_explorer
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(json_file), intent(inout) :: exploring_data

        class(Abstract_Maximum_Box_Compression), allocatable :: maximum_box_compression
        logical :: measure_maximum_box_compression, measure_inv_pow_activities
        logical, dimension(size(physical_model%mixture%components)) :: can_exchange, &
            have_positions, have_orientations

        measure_maximum_box_compression = measure_maximum_compression(exploring_data, &
            volume_change_prefix)
        call maximum_box_compression_create(maximum_box_compression, physical_model%environment%&
            periodic_box, measure_maximum_box_compression)
        call maximum_box_compression_explorer_create(markov_chain_explorer%&
            maximum_box_compression_explorer, physical_model, maximum_box_compression, &
            measure_maximum_box_compression)
        call maximum_box_compression_destroy(maximum_box_compression)
        call volume_change_method_create(markov_chain_explorer%volume_change_method, &
            physical_model, measure_pressure(exploring_data, volume_change_prefix))
        measure_inv_pow_activities = measure_chemical_potentials(exploring_data, &
            particle_insertion_prefix)
        can_exchange = measure_inv_pow_activities ! as if exchange
        call set_have_positions(have_positions, physical_model%mixture%components)
        call box_create(markov_chain_explorer%particle_insertion_domain, physical_model%&
            environment%periodic_box, physical_model%environment%visitable_walls, any(can_exchange)&
            .and. any(have_positions), exploring_data, particle_insertion_prefix)
        call random_coordinates_create(markov_chain_explorer%random_position, &
            markov_chain_explorer%particle_insertion_domain, can_exchange, have_positions)
        call set_have_orientations(have_orientations, physical_model%mixture%components)
        call random_coordinates_create(markov_chain_explorer%random_orientation, can_exchange, &
            have_orientations)
        call particle_insertion_method_create(markov_chain_explorer%particle_insertion_method, &
            physical_model, markov_chain_explorer%random_position, markov_chain_explorer%&
            random_orientation, measure_inv_pow_activities, exploring_data, &
            particle_insertion_prefix)
    end subroutine create

    subroutine destroy(markov_chain_explorer)
        type(Markov_Chain_Explorer_Wrapper), intent(inout) :: markov_chain_explorer

        call particle_insertion_method_destroy(markov_chain_explorer%particle_insertion_method)
        call random_coordinates_destroy(markov_chain_explorer%random_orientation)
        call random_coordinates_destroy(markov_chain_explorer%random_position)
        call box_destroy(markov_chain_explorer%particle_insertion_domain)
        call volume_change_method_destroy(markov_chain_explorer%volume_change_method)
        call maximum_box_compression_explorer_destroy(markov_chain_explorer%&
            maximum_box_compression_explorer)
    end subroutine destroy

end module procedures_markov_chain_explorer_factory
