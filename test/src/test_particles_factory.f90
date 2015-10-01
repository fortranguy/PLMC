program test_particles_factory

use, intrinsic :: iso_fortran_env, only: output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists
use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_factory, only: environment_factory_create, environment_factory_destroy
use types_particles_wrapper, only: Mixture_Wrapper
use procedures_particles_factory, only: particles_factory_create, particles_factory_destroy
implicit none

    type(Environment_Wrapper) :: environment
    type(Mixture_Wrapper) :: mixture

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()
    data_filename = "particles_factory.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call environment_factory_create(environment%periodic_box, input_data, "Environment.")
    call environment_factory_create(environment%floor_penetration, input_data, "Environment.")
    call particles_factory_create(mixture%components(1), input_data, "Mixture.Component 1.", &
        environment)
    call particles_factory_create(mixture%components(2), input_data, "Mixture.Component 2.", &
        environment)
    call particles_factory_create(mixture%inter_diameter, mixture%components(1)%diameter, &
        mixture%components(2)%diameter, input_data, "Mixture.Inter 12.")
    write(*, *) "inter diameter", mixture%inter_diameter%get()
    write(*, *) "minimum inter diameter", mixture%inter_diameter%get_min()

    call particles_factory_destroy(mixture%inter_diameter)
    call particles_factory_destroy(mixture%components(2))
    call particles_factory_destroy(mixture%components(1))
    call environment_factory_destroy(environment%floor_penetration)
    call environment_factory_destroy(environment%periodic_box)
    call input_data%destroy()

end program test_particles_factory
