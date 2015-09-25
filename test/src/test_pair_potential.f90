program test_pair_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use procedures_errors, only: error_exit
use class_particles_diameter, only: Abstract_Particles_Diameter, Concrete_Particles_Diameter
use procedures_particles_factory, only: allocate_and_set_diameter
use class_potential_expression, only: Abstract_Potential_Expression
use class_pair_potential, only: Abstract_Pair_Potential
use procedures_short_potential_factory, only: allocate_and_set_expression, &
    allocate_and_construct_pair

implicit none

    class(Abstract_Pair_Potential), allocatable :: pair_potential
    class(Abstract_Potential_Expression), allocatable :: potential_expression
    class(Abstract_Particles_Diameter), allocatable :: particles_diameter

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found

    real(DP) :: energy, distance
    logical :: overlap

    call json_initialize()
    data_filename = "pair_potential.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call allocate_and_set_diameter(particles_diameter, input_data, "Test Pair Potential.Particles")
    call allocate_and_set_expression(potential_expression, input_data, &
        "Test Pair Potential.Particles", particles_diameter)
    call allocate_and_construct_pair(pair_potential, input_data, "Test Pair Potential.Particles", &
        particles_diameter, potential_expression)

    data_field = "Test Pair Potential.Particles.Potential.distance"
    call input_data%get(data_field, distance, data_found)
    call test_data_found(data_field, data_found)
    call pair_potential%meet(overlap, energy, distance)
    if (overlap) then
        write(output_unit, *) "overlap"
    else
        write(output_unit, *) "energy =", energy
    end if
    write(output_unit, *) "maximum distance", pair_potential%get_max_distance()

    call pair_potential%destroy()
    deallocate(pair_potential)
    deallocate(potential_expression)
    deallocate(particles_diameter)
    call input_data%destroy()

end program test_pair_potential
