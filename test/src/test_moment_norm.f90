module procedures_moment_norm_manipulate

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file
use procedures_checks, only: check_data_found
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm

implicit none

private
public manipulate_moment_norm

contains

    subroutine manipulate_moment_norm(object_name, moment_norm, input_data)
        character(len=*), intent(in) :: object_name
        class(Abstract_Particles_Moment_Norm), intent(inout) :: moment_norm
        type(json_file) :: input_data

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: moment_norm_value

        write(output_unit, *) object_name
        data_field = "Moment Norm.value"
        call input_data%get(data_field, moment_norm_value, data_found)
        call check_data_found(data_field, data_found)
        call moment_norm%set(moment_norm_value)
        write(output_unit, *) "norm =", moment_norm%get()

        deallocate(data_field)
    end subroutine manipulate_moment_norm

end module procedures_moment_norm_manipulate

program test_moment_norm

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use procedures_checks, only: check_file_exists
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm, Null_Particles_Moment_Norm, Concrete_Particles_Moment_Norm
use procedures_moment_norm_manipulate, only: manipulate_moment_norm

implicit none

    class(Abstract_Particles_Moment_Norm), allocatable :: moment_norm
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()

    data_filename = "moment_norm.json"
    call check_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(Null_Particles_Moment_Norm :: moment_norm)
    call manipulate_moment_norm("Null", moment_norm, input_data)
    deallocate(moment_norm)

    allocate(Concrete_Particles_Moment_Norm :: moment_norm)
    call manipulate_moment_norm("Uniform", moment_norm, input_data)
    deallocate(moment_norm)

    deallocate(data_filename)
    call input_data%destroy()

end program test_moment_norm
