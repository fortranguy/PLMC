module procedures_exploring_writers_factory

use data_input_prefixes, only: writers_prefix
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_number_to_string, only: Concrete_Number_to_String
use types_string_wrapper, only: String_Wrapper
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use classes_pair_potential, only: Pair_Potential_Wrapper, Pair_Potentials_Line
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

    subroutine create(writers, environment, wall_pairs, components, short_pairs, &
        volume_change_method, particle_insertion_method, visit_energies, generating_data)
        type(Exploring_Writers_Wrapper), intent(out) :: writers
        type(Environment_Wrapper), intent(in) :: environment
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        type(Pair_Potentials_Line), intent(in) :: short_pairs(:)
        class(Abstract_Volume_Change_Method), intent(in) :: volume_change_method
        class(Abstract_Particle_Insertion_Method), intent(in) :: particle_insertion_method
        logical, intent(in) :: visit_energies
        type(json_file), intent(inout) :: generating_data

        type(String_Wrapper), dimension(size(environment%periodic_boxes)) :: boxes_path
        logical :: selectors(size(components, 1), size(components, 2)), measure_pressure
        character(len=:), allocatable :: make_directory_cmd, separator
        logical :: data_found
        character(len=:), allocatable :: data_field
        type(Concrete_Number_to_String) :: string
        integer :: i_box, box_stat_i

        data_field = writers_prefix//"path separator"
        call generating_data%get(data_field, separator, data_found)
        call check_data_found(data_field, data_found)
        data_field = writers_prefix//"make directory command"
        call generating_data%get(data_field, make_directory_cmd, data_found)
        call check_data_found(data_field, data_found)
        do i_box = 1, size(boxes_path)
            boxes_path(i_box)%string = "exploration_box_"//string%get(i_box)//separator
            call execute_command_line(make_directory_cmd//" "//boxes_path(i_box)%string, &
                exitstat=box_stat_i)
            if (box_stat_i /= 0) call error_exit("procedures_exploring_writers_factory: create:"//&
                " directory can't be created.")
        end do

        measure_pressure = property_measure_pressure(volume_change_method)
        call real_writer_create(writers%maximum_boxes_compression_delta, boxes_path, &
            "maximum_box_compression_delta.out", measure_pressure)
        call real_writer_create(writers%beta_pressures_excess, boxes_path, "beta_pressure_excess.out", &
            measure_pressure)
        call energies_writers_create(writers%gemc_energies, boxes_path, environment%external_fields, wall_pairs, &
            components, short_pairs, visit_energies)
        selectors = measure_chemical_potentials(particle_insertion_method)
        call line_writer_create(writers%gemc_insertion_successes, boxes_path, "insertion_successes.out", selectors)
        call line_writer_create(writers%gemc_inv_pow_activities, boxes_path, "inv_pow_activities.out", selectors)
    end subroutine create

    subroutine destroy(writers)
        type(Exploring_Writers_Wrapper), intent(inout) :: writers

        call energies_writers_destroy(writers%gemc_energies)
        call line_writer_destroy(writers%gemc_inv_pow_activities)
        call line_writer_destroy(writers%gemc_insertion_successes)
        call real_writer_destroy(writers%beta_pressures_excess)
        call real_writer_destroy(writers%maximum_boxes_compression_delta)
    end subroutine destroy

end module procedures_exploring_writers_factory
