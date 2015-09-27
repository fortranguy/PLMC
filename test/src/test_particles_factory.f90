program test_particles_factory

use, intrinsic :: iso_fortran_env, only: output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists
use class_periodic_box, only: Abstract_Periodic_Box
use procedures_box_factory, only: box_factory_create, box_factory_destroy
use types_particles, only: Mixture_Wrapper
use procedures_particles_factory, only: particles_factory_create, particles_factory_destroy
implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    type(Mixture_Wrapper) :: mixture

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()
    data_filename = "particles_factory.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call box_factory_create(periodic_box, input_data, "Box.")
    call particles_factory_create(mixture%components(1), input_data, "Mixture.Component 1.", &
        periodic_box)
    call particles_factory_create(mixture%components(2), input_data, "Mixture.Component 2.", &
        periodic_box)
    call particles_factory_create(mixture%inter_diameters, mixture%components(1)%diameter, &
        mixture%components(2)%diameter, input_data, "Mixture.Inter 12.")
    write(*, *) "inter diameter", mixture%inter_diameters%get()
    write(*, *) "minimum inter diameter", mixture%inter_diameters%get_min()

    call particles_factory_destroy(mixture%inter_diameters)
    call particles_factory_destroy(mixture%components(2))
    call particles_factory_destroy(mixture%components(1))
    call box_factory_destroy(periodic_box)
    call input_data%destroy()

end program test_particles_factory
