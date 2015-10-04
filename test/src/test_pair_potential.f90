program test_pair_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use procedures_checks, only: check_file_exists, check_data_found
use procedures_errors, only: error_exit
use class_particles_diameter, only: Abstract_Particles_Diameter, Concrete_Particles_Diameter
use procedures_particles_factory, only: particles_factory_create, particles_factory_destroy
use class_potential_expression, only: Abstract_Potential_Expression
use class_pair_potential, only: Abstract_Pair_Potential
use procedures_short_potential_factory, only: short_potential_factory_create, &
    short_potential_factory_destroy

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
    call check_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call particles_factory_create(particles_diameter, input_data, "Test Pair Potential.Particles")
    call short_potential_factory_create(potential_expression, input_data, &
        "Test Pair Potential.Particles", particles_diameter)
    call short_potential_factory_create(pair_potential, input_data, &
        "Test Pair Potential.Particles", particles_diameter, potential_expression)

    data_field = "Test Pair Potential.Particles.Potential.distance"
    call input_data%get(data_field, distance, data_found)
    call check_data_found(data_field, data_found)
    call pair_potential%meet(overlap, energy, distance)
    if (overlap) then
        write(output_unit, *) "overlap"
    else
        write(output_unit, *) "energy =", energy
    end if
    write(output_unit, *) "maximum distance", pair_potential%get_max_distance()

    call pair_potential%destroy()
    deallocate(pair_potential)
    call short_potential_factory_destroy(potential_expression)
    call particles_factory_destroy(particles_diameter)
    call input_data%destroy()

end program test_pair_potential
