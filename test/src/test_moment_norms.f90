module procedures_moment_norms_manipulate

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file
use module_data, only: test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_moment_norms, only: Abstract_Moment_Norms

implicit none

private
public manipulate_moment_norms

contains

    subroutine manipulate_moment_norms(object_name, moment_norms, input_data)
    
        character(len=*), intent(in) :: object_name
        class(Abstract_Moment_Norms), intent(inout) :: moment_norms
        type(json_file) :: input_data
        
        class(Abstract_Particles_Number), allocatable :: particles_num
        character(len=:), allocatable :: data_field
        logical :: found
        real(DP) :: moment_norms_norm
        integer :: moment_norms_num, i_particle
    
        allocate(Concrete_Particles_Number :: particles_num)
        data_field = "Moment Norms.number"
        call input_data%get(data_field, moment_norms_num, found)
        call test_data_found(data_field, found)
        call particles_num%set(moment_norms_num)
    
        write(output_unit, *) object_name
        call moment_norms%construct(particles_num)
        data_field = "Moment Norms.norm"
        call input_data%get(data_field, moment_norms_norm, found)
        call test_data_found(data_field, found)
        call moment_norms%set(1, moment_norms_norm)
        do i_particle = 1, particles_num%get()
            write(output_unit, *) "norm", i_particle, "=", moment_norms%get(i_particle)
        end do
        call moment_norms%destroy()
        
        deallocate(data_field)
        deallocate(particles_num)
        
    end subroutine manipulate_moment_norms

end module procedures_moment_norms_manipulate

program test_moment_norms

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists
use class_moment_norms, only: Abstract_Moment_Norms, Null_Moment_Norms, Uniform_Moment_Norms
use procedures_moment_norms_manipulate, only: manipulate_moment_norms

implicit none
    
    class(Abstract_Moment_Norms), allocatable :: moment_norms
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()

    data_filename = "moment_norms.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Null_Moment_Norms :: moment_norms)
    call manipulate_moment_norms("Null", moment_norms, input_data)
    deallocate(moment_norms)
    
    allocate(Uniform_Moment_Norms :: moment_norms)
    call manipulate_moment_norms("Uniform", moment_norms, input_data)
    deallocate(moment_norms)

    deallocate(data_filename)
    call input_data%destroy()

end program test_moment_norms
