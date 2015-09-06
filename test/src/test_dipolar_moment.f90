program test_dipolar_moment

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit

use module_data, only: test_file_exists, test_data_found
use class_dipolar_moment, only : Abstract_Dipolar_Moment, Concrete_Dipolar_Moment
use json_module, only: json_file, json_initialize

implicit none

    class(Abstract_Dipolar_Moment), allocatable :: dipolar_moment
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    real(DP) :: norm

    call json_initialize()

    data_filename = "dipolar_moment.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Concrete_Dipolar_Moment :: dipolar_moment)

    data_field = "Dipolar_Moment.norm"
    call input_data%get(data_field, norm, found)
    call test_data_found(data_field, found)

    call dipolar_moment%set_norm(norm)
    write(output_unit, *) "norm =", dipolar_moment%get_norm()

    deallocate(dipolar_moment)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_dipolar_moment