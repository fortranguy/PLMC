program test_spheres

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit

use module_data, only: test_file_exists, test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_spheres, only: Abstract_Spheres, Uniform_Spheres
use json_module, only: json_file, json_initialize

implicit none
    
    class(Abstract_Particles_Number), allocatable :: particles_num
    class(Abstract_Spheres), allocatable :: spheres
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    integer :: spheres_num, i_particle
    real(DP) :: spheres_diameter

    call json_initialize()

    data_filename = "spheres.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Concrete_Particles_Number :: particles_num)    
    data_field = "Spheres.number"
    call input_data%get(data_field, spheres_num, found)
    call test_data_found(data_field, found)
    call particles_num%set(spheres_num)
    
    allocate(Uniform_Spheres :: spheres)
    call spheres%construct(particles_num)
    data_field = "Spheres.diameter"
    call input_data%get(data_field, spheres_diameter, found)
    call test_data_found(data_field, found)
    call spheres%set_diameter(1, spheres_diameter)
    
    do i_particle = 1, particles_num%get()
        write(output_unit, *) "diameter", i_particle, "=", spheres%get_diameter(i_particle)
    end do

    call spheres%destroy()
    deallocate(particles_num)
    deallocate(spheres)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_spheres
