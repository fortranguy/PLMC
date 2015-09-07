program test_inter_spheres

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit

use module_data, only: test_file_exists, test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_spheres, only: Abstract_Spheres, Uniform_Spheres
use class_inter_spheres, only: Abstract_Inter_Spheres, Concrete_Inter_Spheres
use json_module, only: json_file, json_initialize

implicit none
    
    class(Abstract_Particles_Number), allocatable :: particles_num_1, particles_num_2
    class(Abstract_Spheres), allocatable :: spheres_1, spheres_2
    class(Abstract_Inter_Spheres), allocatable :: inter_spheres_12
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    integer :: spheres_1_num, i_particle
    integer :: spheres_2_num, j_particle
    real(DP) :: spheres_1_diameter, spheres_2_diameter, non_additivity

    call json_initialize()

    data_filename = "inter_spheres.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Concrete_Particles_Number :: particles_num_1)
    data_field = "Spheres 1.number"
    call input_data%get(data_field, spheres_1_num, found)
    call test_data_found(data_field, found)
    call particles_num_1%set_num(spheres_1_num)
    
    allocate(Uniform_Spheres :: spheres_1)
    data_field = "Spheres 1.diameter"
    call input_data%get(data_field, spheres_1_diameter, found)
    call test_data_found(data_field, found)
    call spheres_1%construct(particles_num_1)
    call spheres_1%set_diameter(1, spheres_1_diameter)
    
    allocate(Concrete_Particles_Number :: particles_num_2)
    data_field = "Spheres 2.number"
    call input_data%get(data_field, spheres_2_num, found)
    call test_data_found(data_field, found)
    call particles_num_2%set_num(spheres_2_num)
    
    allocate(Uniform_Spheres :: spheres_2)
    data_field = "Spheres 2.diameter"
    call input_data%get(data_field, spheres_2_diameter, found)
    call test_data_found(data_field, found)
    call spheres_2%construct(particles_num_2)
    call spheres_2%set_diameter(1, spheres_2_diameter)
    
    allocate(Concrete_Inter_Spheres :: inter_spheres_12)    
    data_field = "Inter Spheres 12.non additivity"
    call input_data%get(data_field, non_additivity, found)
    call test_data_found(data_field, found)
    call inter_spheres_12%construct(spheres_1, spheres_2, non_additivity)
    
    do i_particle = 1, particles_num_1%get_num()
        do j_particle = 1, particles_num_2%get_num()
            write(output_unit, *) "diameter", i_particle, j_particle, &
                inter_spheres_12%get_diameter(i_particle, j_particle)
        end do
    end do

    call inter_spheres_12%destroy()
    deallocate(inter_spheres_12)
    call spheres_1%destroy()
    deallocate(spheres_1)
    call spheres_2%destroy()
    deallocate(spheres_2)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_inter_spheres
