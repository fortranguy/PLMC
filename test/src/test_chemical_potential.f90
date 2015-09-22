program test_chemical_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_particles_chemical_potential, only: Abstract_Particles_Chemical_Potential, &
    Concrete_Particles_Chemical_Potential

implicit none

    class(Abstract_Particles_Chemical_Potential), allocatable :: chemical_potential
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    real(DP) :: density, excess

    call json_initialize()
    data_filename = "chemical_potential.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(Concrete_Particles_Chemical_Potential :: chemical_potential)
    data_field = "Chemical Potential.density"
    call input_data%get(data_field, density, found)
    call test_data_found(data_field, found)
    data_field = "Chemical Potential.excess"
    call input_data%get(data_field, excess, found)
    call test_data_found(data_field, found)
    call chemical_potential%set(density, excess)

    write(output_unit, *) "density =", chemical_potential%get_density()
    write(output_unit, *) "excess =", chemical_potential%get_excess()

    deallocate(chemical_potential)
    call input_data%destroy()

end program test_chemical_potential
