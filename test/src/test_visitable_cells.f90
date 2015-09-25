program test_visitable_cells

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use procedures_box_factory, only: allocate_and_set_periodic_box
use class_particles_number, only: Abstract_Particles_Number
use class_particles_diameter, only: Abstract_Particles_Diameter
use class_particles_positions, only: Abstract_Particles_Positions
use procedures_particles_factory, only: allocate_and_set_number, allocate_and_set_diameter, &
    allocate_and_construct_positions, set_positions
use class_potential_expression, only: Abstract_Potential_Expression
use class_pair_potential, only: Abstract_Pair_Potential
use class_particles_potential, only: Abstract_Particles_Potential
use procedures_short_potential_factory, only: allocate_and_set_expression, &
    allocate_and_construct_pair, allocate_and_construct_particles_potential
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
    character(len=:), allocatable :: data_filename, data_field, list_data_structure
    logical :: data_found

    type(Concrete_Particle) :: particle
    real(DP) :: energy, energy_i
    integer :: i_particle
    logical :: overlap

    call json_initialize()
    data_filename = "particles_potential.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call allocate_and_set_periodic_box(periodic_box, input_data, "Test Particles Potential")
    call allocate_and_set_number(particles_number, input_data, &
        "Test Particles Potential.Particles")
    call allocate_and_set_diameter(particles_diameter, input_data, &
        "Test Particles Potential.Particles")
    call allocate_and_construct_positions(particles_positions, periodic_box, particles_number)
    call set_positions(particles_positions, input_data, "Test Particles Potential.Particles")
    call allocate_and_set_expression(potential_expression, input_data, &
        "Test Particles Potential.Particles", particles_diameter)
    call allocate_and_construct_pair(pair_potential, input_data, &
        "Test Particles Potential.Particles", particles_diameter, potential_expression)

    data_field = "Test Particles Potential.Memory.data structure"
    call input_data%get(data_field, list_data_structure, data_found)
    call test_data_found(data_field, data_found)
    select case(list_data_structure)
        case ("list")
            allocate(Concrete_Visitable_List :: visitable_list)
        case ("array")
            allocate(Concrete_Visitable_Array :: visitable_list)
        case default
            call error_exit(list_data_structure//" unknown.")
    end select
    deallocate(list_data_structure)

    select type(periodic_box)
        type is (XYZ_Periodic_Box)
            allocate(XYZ_PBC_Visitable_Cells :: visitable_cells)
        type is (XY_Periodic_Box)
            allocate(XY_PBC_Visitable_Cells :: visitable_cells)
        class default
            call error_exit("Periodic_Box type unknown.")
    end select

    call visitable_cells%construct(visitable_list, periodic_box, particles_positions, &
        pair_potential)
    call visitable_cells%fill()

    energy = 0._DP
    do i_particle = 1, particles_positions%get_num()
        particle%same_type = .true.
        particle%i = i_particle
        particle%position = particles_positions%get(particle%i)
        call visitable_cells%visit(particle, overlap, energy_i)
        if (overlap) exit
        energy = energy + energy_i
    end do
    if (overlap) then
        write(output_unit,*) "overlap"
    else
        energy = energy / 2._DP
        write(output_unit, *) "energy =", energy
    end if

    call visitable_cells%destroy()
    deallocate(visitable_cells)
    deallocate(visitable_list)
    call pair_potential%destroy()
    deallocate(pair_potential)
    deallocate(potential_expression)
    call particles_positions%destroy()
    deallocate(particles_positions)
    deallocate(particles_diameter)
    deallocate(particles_number)
    deallocate(periodic_box)
    call input_data%destroy()

end program test_visitable_cells
