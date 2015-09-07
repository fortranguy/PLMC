program test_dipolar_moment_norms

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit

use module_data, only: test_file_exists, test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_dipolar_moment_norms, only: Abstract_Dipolar_Moment_Norms, Uniform_Dipolar_Moment_Norms
use json_module, only: json_file, json_initialize

implicit none
    
    class(Abstract_Particles_Number), allocatable :: particles_num
    class(Abstract_Dipolar_Moment_Norms), allocatable :: dipolar_moment_norms
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    integer :: dipolar_moment_norms_num, i_particle
    real(DP) :: dipolar_moment_norms_norm

    call json_initialize()

    data_filename = "dipolar_moment_norms.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Concrete_Particles_Number :: particles_num)    
    data_field = "Dipolar Moment Norms.number"
    call input_data%get(data_field, dipolar_moment_norms_num, found)
    call test_data_found(data_field, found)
    call particles_num%set_num(dipolar_moment_norms_num)
    
    allocate(Uniform_Dipolar_Moment_Norms :: dipolar_moment_norms)
    call dipolar_moment_norms%construct(particles_num)
    data_field = "Dipolar Moment Norms.norm"
    call input_data%get(data_field, dipolar_moment_norms_norm, found)
    call test_data_found(data_field, found)
    call dipolar_moment_norms%set_norm(1, dipolar_moment_norms_norm)
    
    do i_particle = 1, particles_num%get_num()
        write(output_unit, *) "norm", i_particle, "=", dipolar_moment_norms%get_norm(i_particle)
    end do

    call dipolar_moment_norms%destroy()
    deallocate(particles_num)
    deallocate(dipolar_moment_norms)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_dipolar_moment_norms
