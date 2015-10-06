program test_particles_factory

use, intrinsic :: iso_fortran_env, only: output_unit
use json_module, only: json_file, json_initialize
use procedures_checks, only: check_file_exists
use procedures_property_inquirers, only: particles_exist
use class_periodic_box, only: Abstract_Periodic_Box
use procedures_environment_factory, only: environment_factory_create, environment_factory_destroy
use types_particles_wrapper, only: Mixture_Wrapper
use procedures_particles_factory, only: particles_factory_create, particles_factory_destroy
implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    type(Mixture_Wrapper) :: mixture

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename
    logical :: mixture_exists

    call json_initialize()
    data_filename = "particles_factory.json"
    call check_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call environment_factory_create(periodic_box, input_data, "Environment.")
    call particles_factory_create(mixture%components(1), periodic_box, input_data, &
        "Mixture.Component 1.")
    call particles_factory_create(mixture%components(2), periodic_box, input_data, &
        "Mixture.Component 2.")
    mixture_exists = particles_exist(mixture%components(1)%number) .and. &
            particles_exist(mixture%components(2)%number)
    call particles_factory_create(mixture%inter_diameter, mixture_exists, &
        mixture%components(1)%diameter, mixture%components(2)%diameter, input_data, &
        "Mixture.Inter 12.")
    write(*, *) "inter diameter", mixture%inter_diameter%get()
    write(*, *) "minimum inter diameter", mixture%inter_diameter%get_min()

    call particles_factory_destroy(mixture%inter_diameter)
    call particles_factory_destroy(mixture%components(2))
    call particles_factory_destroy(mixture%components(1))
    call environment_factory_destroy(periodic_box)
    call input_data%destroy()

end program test_particles_factory
