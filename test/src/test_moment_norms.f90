program test_moment_norms

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_moment_norms, only: Abstract_Moment_Norms, Uniform_Moment_Norms

implicit none
    
    class(Abstract_Particles_Number), allocatable :: particles_num
    class(Abstract_Moment_Norms), allocatable :: moment_norms
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    integer :: moment_norms_num, i_particle
    real(DP) :: moment_norms_norm

    call json_initialize()

    data_filename = "moment_norms.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Concrete_Particles_Number :: particles_num)
    data_field = "Moment Norms.number"
    call input_data%get(data_field, moment_norms_num, found)
    call test_data_found(data_field, found)
    call particles_num%set(moment_norms_num)
    
    allocate(Uniform_Moment_Norms :: moment_norms)
    call moment_norms%construct(particles_num)
    data_field = "Moment Norms.norm"
    call input_data%get(data_field, moment_norms_norm, found)
    call test_data_found(data_field, found)
    call moment_norms%set(1, moment_norms_norm)
    
    do i_particle = 1, particles_num%get()
        write(output_unit, *) "norm", i_particle, "=", moment_norms%get(i_particle)
    end do

    call moment_norms%destroy()
    deallocate(particles_num)
    deallocate(moment_norms)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_moment_norms
