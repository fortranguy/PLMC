program test_environment_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists
use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_factory, only: environment_factory_create, environment_factory_destroy

implicit none

    type(Environment_Wrapper) :: environment

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()
    data_filename = "environment_factory.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call environment_factory_create(environment, input_data, "Environment.")
    write(output_unit, *) "box size", environment%periodic_box%get_size()
    write(output_unit, *) "temperature", environment%temperature%get()
    write(output_unit, *) "external field at origin", &
        environment%external_field%get([0._DP, 0._DP, 0._DP])
    write(output_unit, *) "reciprocal lattice numbers", environment%reciprocal_lattice%get_numbers()
    write(output_unit, *) "walls gap", environment%walls_potential%get_gap()

    call environment_factory_destroy(environment)
    call input_data%destroy()

end program test_environment_factory
