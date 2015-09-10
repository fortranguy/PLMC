module procedures_moments_norm_manipulate

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file
use module_data, only: test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_moments_norm, only: Abstract_Moments_Norm

implicit none

private
public manipulate_moments_norm

contains

    subroutine manipulate_moments_norm(object_name, moments_norm, input_data)
        character(len=*), intent(in) :: object_name
        class(Abstract_Moments_Norm), intent(inout) :: moments_norm
        type(json_file) :: input_data
        
        class(Abstract_Particles_Number), allocatable :: particles_number
        character(len=:), allocatable :: data_field
        logical :: found
        real(DP) :: moments_norm_value
        integer :: moments_norm_num, i_particle
    
        allocate(Concrete_Particles_Number :: particles_number)
        data_field = "Moments Norm.number"
        call input_data%get(data_field, moments_norm_num, found)
        call test_data_found(data_field, found)
        call particles_number%set(moments_norm_num)
    
        write(output_unit, *) object_name
        call moments_norm%construct(particles_number)
        data_field = "Moments Norm.value"
        call input_data%get(data_field, moments_norm_value, found)
        call test_data_found(data_field, found)
        call moments_norm%set(1, moments_norm_value)
        do i_particle = 1, particles_number%get()
            write(output_unit, *) "norm", i_particle, "=", moments_norm%get(i_particle)
        end do
        data_field = "Moments Norm.add"
        call input_data%get(data_field, moments_norm_value, found)
        call test_data_found(data_field, found)
        call particles_number%set(particles_number%get() + 1)
        call moments_norm%add(moments_norm_value)
        write(output_unit, *) "new norm", particles_number%get(), "=", &
            moments_norm%get(particles_number%get())
        
        data_field = "Moments Norm.remove"
        call input_data%get(data_field, i_particle, found)
        call test_data_found(data_field, found)
        call moments_norm%remove(i_particle)
        call particles_number%set(particles_number%get() - 1)
        do i_particle = 1, particles_number%get()
            write(output_unit, *) "norm", i_particle, "=", moments_norm%get(i_particle)
        end do
        call moments_norm%destroy()
        
        deallocate(data_field)
        deallocate(particles_number)
        
    end subroutine manipulate_moments_norm

end module procedures_moments_norm_manipulate

program test_moments_norm

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists
use class_moments_norm, only: Abstract_Moments_Norm, Null_Moments_Norm, Uniform_Moments_Norm
use procedures_moments_norm_manipulate, only: manipulate_moments_norm

implicit none
    
    class(Abstract_Moments_Norm), allocatable :: moments_norm
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()

    data_filename = "moments_norm.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Null_Moments_Norm :: moments_norm)
    call manipulate_moments_norm("Null", moments_norm, input_data)
    deallocate(moments_norm)
    
    allocate(Uniform_Moments_Norm :: moments_norm)
    call manipulate_moments_norm("Uniform", moments_norm, input_data)
    deallocate(moments_norm)

    deallocate(data_filename)
    call input_data%destroy()

end program test_moments_norm
