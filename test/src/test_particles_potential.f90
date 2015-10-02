program test_particles_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box
use procedures_environment_factory, only: environment_factory_create, environment_factory_destroy
use class_particles_number, only: Abstract_Particles_Number
use class_particles_diameter, only: Abstract_Particles_Diameter
use class_particles_positions, only: Abstract_Particles_Positions
use procedures_particles_factory, only: particles_factory_create, particles_factory_set, &
    particles_factory_destroy
use class_potential_expression, only: Abstract_Potential_Expression
use class_pair_potential, only: Abstract_Pair_Potential
use class_particles_potential, only: Abstract_Particles_Potential
use procedures_short_potential_factory, only: short_potential_factory_create, &
    short_potential_factory_destroy
use types_particle, only: Concrete_Particle
use class_particles_potential, only: Abstract_Particles_Potential

implicit none

    class(Abstract_Particles_Potential), allocatable :: particles_potential
    type(Concrete_Particle) :: particle
    class(Abstract_Pair_Potential), allocatable :: pair_potential
    class(Abstract_Potential_Expression), allocatable :: potential_expression
    class(Abstract_Particles_Positions), allocatable :: particles_positions
    class(Abstract_Particles_Diameter), allocatable :: particles_diameter
    class(Abstract_Particles_Number), allocatable :: particles_number
    class(Abstract_Periodic_Box), allocatable :: periodic_box

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    real(DP) :: energy, energy_i
    integer :: i_particle
    logical :: overlap

    call json_initialize()
    data_filename = "particles_potential.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call environment_factory_create(periodic_box, input_data, "Box.")
    call particles_factory_create(particles_number, input_data, "Particles.")
    call particles_factory_create(particles_diameter, input_data, "Particles.")
    call particles_factory_create(particles_positions, periodic_box, particles_number)
    call particles_factory_set(particles_positions, input_data, "Particles.")
    call short_potential_factory_create(potential_expression, input_data, &
        "Particles.Potential.", particles_diameter)
    call short_potential_factory_create(pair_potential, input_data, &
        "Particles.Potential.", particles_diameter, potential_expression)
    call short_potential_factory_create(particles_potential, periodic_box, particles_positions)

    energy = 0._DP
    do i_particle = 1, particles_positions%get_num()
        particle%same_type = .true.
        particle%i = i_particle
        particle%position = particles_positions%get(particle%i)
        call particles_potential%visit(overlap, energy_i, particle, pair_potential)
        if (overlap) exit
        energy = energy + energy_i
    end do
    if (overlap) then
        write(output_unit,*) "overlap"
    else
        energy = energy / 2._DP
        write(output_unit, *) "energy =", energy
    end if

    call particles_potential%destroy()
    deallocate(particles_potential)
    call pair_potential%destroy()
    deallocate(pair_potential)
    call short_potential_factory_destroy(potential_expression)
    call particles_factory_destroy(particles_positions)
    call particles_factory_destroy(particles_diameter)
    call particles_factory_destroy(particles_number)
    call environment_factory_destroy(periodic_box)
    call input_data%destroy()

end program test_particles_potential
