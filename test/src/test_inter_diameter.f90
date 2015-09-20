program test_inter_diameter

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_particles_diameter, only: Abstract_Particles_Diameter, Concrete_Particles_Diameter
use class_inter_particles_diameter, only: Abstract_Inter_Particles_Diameter, &
    Concrete_Inter_Particles_Diameter

implicit none

    class(Abstract_Particles_Diameter), allocatable :: diameter_1, diameter_2
    class(Abstract_Inter_Particles_Diameter), allocatable :: inter_diameter_12
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found
    real(DP) :: diameter_1_value, diameter_1_min_factor
    real(DP) :: diameter_2_value, diameter_2_min_factor
    real(DP) :: inter_diameter_12_offset

    call json_initialize()

    data_filename = "inter_diameter.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    allocate(Concrete_Particles_Diameter :: diameter_1)
    data_field = "Diameter 1.value"
    call input_data%get(data_field, diameter_1_value, data_found)
    call test_data_found(data_field, data_found)
    data_field = "Diameter 1.minimum factor"
    call input_data%get(data_field, diameter_1_min_factor, data_found)
    call test_data_found(data_field, data_found)
    call diameter_1%set(diameter_1_value, diameter_1_min_factor)

    allocate(Concrete_Particles_Diameter :: diameter_2)
    data_field = "Diameter 2.value"
    call input_data%get(data_field, diameter_2_value, data_found)
    call test_data_found(data_field, data_found)
    data_field = "Diameter 2.minimum factor"
    call input_data%get(data_field, diameter_2_min_factor, data_found)
    call test_data_found(data_field, data_found)
    call diameter_2%set(diameter_2_value, diameter_2_min_factor)

    allocate(Concrete_Inter_Particles_Diameter :: inter_diameter_12)
    data_field = "Inter Diameter 12.offset"
    call input_data%get(data_field, inter_diameter_12_offset, data_found)
    call test_data_found(data_field, data_found)
    call inter_diameter_12%construct(diameter_1, diameter_2, inter_diameter_12_offset)
    deallocate(data_field)

    write(output_unit, *) "diameter", inter_diameter_12%get()
    write(output_unit, *) "minimum diameter", inter_diameter_12%get_min()

    call inter_diameter_12%destroy()
    deallocate(inter_diameter_12)
    deallocate(diameter_2)
    deallocate(diameter_1)
    call input_data%destroy()

end program test_inter_diameter
