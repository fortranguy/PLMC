program test_inter_diameters

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_diameters, only: Abstract_Diameters, Uniform_Diameters
use class_inter_diameters, only: Abstract_Inter_Diameters, Concrete_Inter_Diameters

implicit none
    
    class(Abstract_Particles_Number), allocatable :: particles_number_1, particles_number_2
    class(Abstract_Diameters), allocatable :: diameters_1, diameters_2
    class(Abstract_Inter_Diameters), allocatable :: inter_diameters_12
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    integer :: diameters_1_num, i_particle
    integer :: diameters_2_num, j_particle
    real(DP) :: diameters_1_diameter, diameters_2_diameter, non_additivity

    call json_initialize()

    data_filename = "inter_diameters.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Concrete_Particles_Number :: particles_number_1)
    data_field = "Diameters 1.number"
    call input_data%get(data_field, diameters_1_num, found)
    call test_data_found(data_field, found)
    call particles_number_1%set(diameters_1_num)
    
    allocate(Uniform_Diameters :: diameters_1)
    data_field = "Diameters 1.diameter"
    call input_data%get(data_field, diameters_1_diameter, found)
    call test_data_found(data_field, found)
    call diameters_1%construct(particles_number_1)
    call diameters_1%set(1, diameters_1_diameter)
    
    allocate(Concrete_Particles_Number :: particles_number_2)
    data_field = "Diameters 2.number"
    call input_data%get(data_field, diameters_2_num, found)
    call test_data_found(data_field, found)
    call particles_number_2%set(diameters_2_num)
    
    allocate(Uniform_Diameters :: diameters_2)
    data_field = "Diameters 2.diameter"
    call input_data%get(data_field, diameters_2_diameter, found)
    call test_data_found(data_field, found)
    call diameters_2%construct(particles_number_2)
    call diameters_2%set(1, diameters_2_diameter)
    
    allocate(Concrete_Inter_Diameters :: inter_diameters_12)
    data_field = "Inter Diameters 12.non additivity"
    call input_data%get(data_field, non_additivity, found)
    call test_data_found(data_field, found)
    call inter_diameters_12%construct(diameters_1, diameters_2, non_additivity)
    
    do i_particle = 1, particles_number_1%get()
        do j_particle = 1, particles_number_2%get()
            write(output_unit, *) "diameter", i_particle, j_particle, &
                inter_diameters_12%get(i_particle, j_particle)
        end do
    end do

    call inter_diameters_12%destroy()
    deallocate(inter_diameters_12)
    call diameters_1%destroy()
    deallocate(diameters_1)
    call diameters_2%destroy()
    deallocate(diameters_2)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_inter_diameters
