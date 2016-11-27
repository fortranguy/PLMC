program test_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use procedures_checks, only: check_file_exists
use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box
use procedures_environment_factory, only: environment_factory_create, environment_factory_destroy
use types_particles_wrapper, only: Particles_Wrapper
use procedures_particles_factory, only: particles_factory_create, particles_factory_destroy
use types_changes_wrapper, only: Changes_Wrapper
use procedures_changes_factory, only: changes_factory_create, changes_factory_destroy

implicit none

    type(Changes_Wrapper) :: changes
    type(Particles_Wrapper) :: particles
    class(Abstract_Periodic_Box), allocatable :: periodic_box

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()
    data_filename = "changes_factory.json"
    call check_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call environment_factory_create(periodic_box, input_data, "Environment.")
    call particles_factory_create(particles, periodic_box, input_data, "Particles.")
    call changes_factory_create(changes, periodic_box, particles, input_data, "Changes.")
    write(output_unit, *) "moved_positions(1)", changes%moved_positions%get(1)
    write(output_unit, *) "rotated_orientations(1)", changes%rotated_orientations%get(1)
    call changes%particles_exchange%remove(particles%number%get())
    write(output_unit, *) "last particle removed"
    write(output_unit, *) "number of particle", particles%number%get()
    write(output_unit, *) "last position", particles%positions%get(particles%number%get())
    write(output_unit, *) "last orientation", particles%orientations%get(particles%number%get())

    call changes_factory_destroy(changes)
    call environment_factory_destroy(periodic_box)
    call input_data%destroy()

end program test_changes_factory
