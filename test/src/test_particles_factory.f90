program test_particles_factory

use, intrinsic :: iso_fortran_env, only: output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists
use class_periodic_box, only: Abstract_Periodic_Box
use procedures_box_factory, only: allocate_and_set_periodic_box
use types_particles, only: Mixture_Wrapper
use procedures_particles_factory, only: mixture_factory_create, mixture_factory_destroy
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

    call allocate_and_set_periodic_box(periodic_box, input_data, "Box")
    call mixture_factory_create(mixture, input_data, "Mixture", periodic_box)
    write(*, *) "inter diameter", mixture%inter_diameters%get()
    write(*, *) "minimum inter diameter", mixture%inter_diameters%get_min()

    call mixture_factory_destroy(mixture)
    deallocate(periodic_box)
    call input_data%destroy()

end program test_particles_factory
