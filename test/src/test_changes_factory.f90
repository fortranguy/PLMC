program test_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use types_particles, only: Particles_Wrapper
use procedures_particles_factory, only: particles_factory_create, particles_factory_destroy
use types_changes, only: Changes_Wrapper
use procedures_changes_factory, only: changes_factory_create, changes_factory_destroy

implicit none

    type(Changes_Wrapper) :: changes
    type(Particles_Wrapper) :: particles
    class(Abstract_Periodic_Box), allocatable :: periodic_box

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field, box_name
    logical :: data_found
    real(DP), allocatable :: box_size(:)

    call json_initialize()
    data_filename = "changes_factory.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    data_field = "Periodic Box.name"
    call input_data%get(data_field, box_name, data_found)
    call test_data_found(data_field, data_found)
    select case(box_name)
        case("XYZ")
            allocate(XYZ_Periodic_Box :: periodic_box)
        case("XY")
            allocate(XY_Periodic_Box :: periodic_box)
        case default
            call error_exit(data_field//" unkown.")
    end select
    deallocate(box_name)
    data_field = "Periodic Box.size"
    call input_data%get(data_field, box_size, data_found)
    call test_data_found(data_field, data_found)
    call periodic_box%set(box_size)

    call particles_factory_create(particles, input_data, "Particles", periodic_box)
    call changes_factory_create(changes, input_data, "Particles", particles)
    write(output_unit, *) "moved_positions(1)", changes%moved_positions%get(1)
    write(output_unit, *) "rotated_orientations(1)", changes%rotated_orientations%get(1)
    call changes%particles_exchange%remove(particles%number%get())
    write(output_unit, *) "last particle removed"
    write(output_unit, *) "number of particle", particles%number%get()
    write(output_unit, *) "last position", particles%positions%get(particles%number%get())
    write(output_unit, *) "last orientation", particles%orientations%get(particles%number%get())

    call changes_factory_destroy(changes)
    call input_data%destroy()
    deallocate(periodic_box)

end program test_changes_factory
