module procedures_markov_chain_explorer_factory

use data_input_prefixes, only: particle_insertion_prefix, volume_change_prefix
use json_module, only: json_file
use procedures_errors, only: warning_continue
use procedures_boxes_factory, only: boxes_create => create, boxes_destroy => destroy
use procedures_mixture_factory, only: set_have_positions, set_have_orientations
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use procedures_random_coordinates_factory, only: random_coordinates_create => create, &
    random_coordinates_destroy => destroy
use classes_changed_box_size_ratio, only: Abstract_Changed_Box_Size_Ratio
use procedures_changed_boxes_size_ratio_factory, only: changed_boxes_size_ratio_create => create, &
    changed_boxes_size_ratio_destroy => destroy
use classes_maximum_box_compression, only: Abstract_Maximum_Box_Compression
use procedures_maximum_box_compression_factory, only: maximum_box_compression_create => create, &
    maximum_box_compression_destroy => destroy
use classes_random_coordinates, only: Abstract_Random_Coordinates
use procedures_maximum_box_compression_explorer_factory, only: &
    maximum_box_compression_explorer_create => create, &
    maximum_box_compression_explorer_destroy => destroy
use procedures_volume_change_method_factory, only: volume_change_method_create => create, &
    volume_change_method_destroy => destroy
use procedures_particle_insertion_method_factory, only: particle_insertion_method_create => create,&
    particle_insertion_method_destroy => destroy
use types_markov_chain_explorer_wrapper, only: Markov_Chain_Explorer_Wrapper
use procedures_exploration_inquirers, only: measure_pressure, measure_chemical_potentials

implicit none

private
public :: create, destroy

contains

    !> @todo
    !> multiple boxes generalisation
    subroutine create(markov_chain_explorer, physical_model, visit_energies, exploring_data)
        type(Markov_Chain_Explorer_Wrapper), intent(out) :: markov_chain_explorer
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        logical, intent(in) :: visit_energies
        type(json_file), intent(inout) :: exploring_data

        class(Abstract_Maximum_Box_Compression), allocatable :: maximum_box_compression
        logical :: boxes_size_can_change(size(physical_model%environment%periodic_boxes))
        logical :: measure_pressure_excess, measure_inv_pow_activities
        logical :: have_positions(size(physical_model%mixture%gemc_components, 1), &
            size(physical_model%mixture%gemc_components, 2))
        logical :: have_orientations(size(have_positions, 1), size(have_positions, 2))
        class(Abstract_Random_Coordinates), allocatable :: random_positions(:), random_orientation
        logical :: can_exchange(size(have_positions, 1), size(have_positions, 2))

        measure_pressure_excess = measure_pressure(exploring_data, volume_change_prefix)
        if (.not.visit_energies .and. measure_pressure_excess) then
            call warning_continue("procedures_markov_chain_explorer_factory: create: "//&
                "measure_pressure_excess needs visit_energies.")
        end if
        call maximum_box_compression_create(maximum_box_compression, physical_model%environment%&
            periodic_boxes, measure_pressure_excess)
        call maximum_box_compression_explorer_create(markov_chain_explorer%&
            maximum_box_compression_explorer, physical_model, maximum_box_compression, &
            measure_pressure_excess)
        call maximum_box_compression_destroy(maximum_box_compression)
        boxes_size_can_change = measure_pressure_excess
        call changed_boxes_size_ratio_create(markov_chain_explorer%changed_boxes_size_ratio, &
            physical_model%environment%periodic_boxes, boxes_size_can_change)
        call volume_change_method_create(markov_chain_explorer%volume_change_method, &
            physical_model, markov_chain_explorer%changed_boxes_size_ratio, measure_pressure_excess, &
            exploring_data)
        measure_inv_pow_activities = measure_chemical_potentials(exploring_data, &
            particle_insertion_prefix)

       call set_have_positions(have_positions, physical_model%mixture%gemc_components)
        can_exchange = measure_inv_pow_activities ! as if exchange

        call boxes_create(markov_chain_explorer%particle_insertion_domains, physical_model%&
            environment%periodic_boxes, physical_model%environment%gemc_visitable_walls, &
            any(can_exchange) .and. any(have_positions), exploring_data, particle_insertion_prefix)
        call random_coordinates_create(random_positions, markov_chain_explorer%&
            particle_insertion_domains, have_positions, can_exchange)
        call set_have_orientations(have_orientations, physical_model%mixture%gemc_components)
        call random_coordinates_create(random_orientation, have_orientations,can_exchange)
        call particle_insertion_method_create(markov_chain_explorer%particle_insertion_method, &
            physical_model, random_positions, random_orientation, measure_inv_pow_activities, &
            exploring_data)
        call random_coordinates_destroy(random_orientation)
        call random_coordinates_destroy(random_positions)
    end subroutine create

    subroutine destroy(markov_chain_explorer)
        type(Markov_Chain_Explorer_Wrapper), intent(inout) :: markov_chain_explorer

        call particle_insertion_method_destroy(markov_chain_explorer%particle_insertion_method)
        call boxes_destroy(markov_chain_explorer%particle_insertion_domains)
        call volume_change_method_destroy(markov_chain_explorer%volume_change_method)
        call changed_boxes_size_ratio_destroy(markov_chain_explorer%changed_boxes_size_ratio)
        call maximum_box_compression_explorer_destroy(markov_chain_explorer%&
            maximum_box_compression_explorer)
    end subroutine destroy

end module procedures_markov_chain_explorer_factory
