program test_visitable_cells

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use procedures_box_factory, only: box_factory_create, box_factory_destroy
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
use class_visitable_list, only: Abstract_Visitable_List, Concrete_Visitable_List, &
    Concrete_Visitable_Array
use class_visitable_cells, only: Abstract_Visitable_Cells, &
    XYZ_PBC_Visitable_Cells, XY_PBC_Visitable_Cells

implicit none

    class(Abstract_Visitable_Cells), allocatable :: visitable_cells
    class(Abstract_Visitable_List), allocatable :: visitable_list
    class(Abstract_Pair_Potential), allocatable :: pair_potential
    class(Abstract_Potential_Expression), allocatable :: potential_expression
    class(Abstract_Particles_Positions), allocatable :: particles_positions
    class(Abstract_Particles_Diameter), allocatable :: particles_diameter
    class(Abstract_Particles_Number), allocatable :: particles_number
    class(Abstract_Periodic_Box), allocatable :: periodic_box

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    type(Concrete_Particle) :: particle
    real(DP) :: energy, energy_i
    integer :: i_particle
    logical :: overlap

    call json_initialize()
    data_filename = "particles_potential.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call box_factory_create(periodic_box, input_data, "Test Particles Potential")
    call particles_factory_create(particles_number, input_data, &
        "Test Particles Potential.Particles")
    call particles_factory_create(particles_diameter, input_data, &
        "Test Particles Potential.Particles")
    call particles_factory_create(particles_positions, periodic_box, particles_number)
    call particles_factory_set(particles_positions, input_data, "Test Particles Potential.Particles")
    call short_potential_factory_create(potential_expression, input_data, &
        "Test Particles Potential.Particles", particles_diameter)
    call short_potential_factory_create(pair_potential, input_data, &
        "Test Particles Potential.Particles", particles_diameter, potential_expression)
    call short_potential_factory_create(visitable_list, input_data, "Test Particles Potential.Particles", &
        particles_positions)
    call short_potential_factory_create(visitable_cells, visitable_list, periodic_box, &
        particles_positions, pair_potential)

    energy = 0._DP
    do i_particle = 1, particles_positions%get_num()
        particle%same_type = .true.
        particle%i = i_particle
        particle%position = particles_positions%get(particle%i)
        call visitable_cells%visit(overlap, energy_i, particle)
        if (overlap) exit
        energy = energy + energy_i
    end do
    if (overlap) then
        write(output_unit,*) "overlap"
    else
        energy = energy / 2._DP
        write(output_unit, *) "energy =", energy
    end if

    call short_potential_factory_destroy(visitable_cells)
    call short_potential_factory_destroy(visitable_list)
    call short_potential_factory_destroy(pair_potential)
    call short_potential_factory_destroy(potential_expression)
    call particles_factory_destroy(particles_positions)
    call particles_factory_destroy(particles_diameter)
    call particles_factory_destroy(particles_number)
    call box_factory_destroy(periodic_box)
    call input_data%destroy()

end program test_visitable_cells
