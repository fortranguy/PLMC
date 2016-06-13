module procedures_markov_chain_explorer_factory

use data_prefixes, only: widom_prefix
use json_module, only: json_file
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use procedures_box_factory, only: box_create => create, box_destroy => destroy
use procedures_mixture_factory, only: set_have_positions, set_have_orientations
use procedures_random_coordinates_factory, only: random_coordinates_create => create, &
    random_coordinates_destroy => destroy
use procedures_widom_method_factory, only: widom_method_create => create, widom_method_destroy => &
    destroy
use types_markov_chain_explorer_wrapper, only: Markov_Chain_Explorer_Wrapper
use procedures_property_inquirers, only: measure_chemical_potentials

implicit none

private
public :: create, destroy

contains

    subroutine create(markov_chain_explorer, physical_model, exploring_data)
        type(Markov_Chain_Explorer_Wrapper), intent(out) :: markov_chain_explorer
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(json_file), intent(inout) :: exploring_data

        logical :: measure_inv_pow_activities
        logical, dimension(size(physical_model%mixture%components)) :: can_exchange, &
            have_positions, have_orientations

        measure_inv_pow_activities = measure_chemical_potentials(exploring_data, widom_prefix)
        can_exchange = measure_inv_pow_activities ! as if exchange
        call set_have_positions(have_positions, physical_model%mixture%components)
        call box_create(markov_chain_explorer%widom_domain, physical_model%environment%&
            periodic_box, physical_model%environment%walls, any(can_exchange) .and. &
            any(have_positions), exploring_data, widom_prefix)
        call random_coordinates_create(markov_chain_explorer%random_position, &
            markov_chain_explorer%widom_domain, can_exchange, have_positions)
        call set_have_orientations(have_orientations, physical_model%mixture%components)
        call random_coordinates_create(markov_chain_explorer%random_orientation, can_exchange, &
            have_orientations)
        call widom_method_create(markov_chain_explorer%widom_method, physical_model, &
            markov_chain_explorer%random_position, markov_chain_explorer%random_orientation, &
            size(can_exchange), measure_inv_pow_activities, exploring_data, widom_prefix)
    end subroutine create

    subroutine destroy(markov_chain_explorer)
        type(Markov_Chain_Explorer_Wrapper), intent(inout) :: markov_chain_explorer

        call widom_method_destroy(markov_chain_explorer%widom_method)
        call random_coordinates_destroy(markov_chain_explorer%random_orientation)
        call random_coordinates_destroy(markov_chain_explorer%random_position)
        call box_destroy(markov_chain_explorer%widom_domain)
    end subroutine destroy

end module procedures_markov_chain_explorer_factory
