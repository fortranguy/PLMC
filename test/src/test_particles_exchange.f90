module types_json_wrapper

use json_module, only: json_value

implicit none

private

    type, public :: JSON_Value_Pointer
        type(json_value), pointer :: ptr
    end type JSON_Value_Pointer

end module types_json_wrapper

module procedures_particles_exchange_write

use json_module, only: json_create_object, json_add, json_print, json_destroy
use types_json_wrapper, only: JSON_Value_Pointer
use types_particles, only: Particles_Wrapper

implicit none

private
public json_write_particles

contains

    subroutine json_write_particles(particles, output_filename)
        type(Particles_Wrapper), intent(in) :: particles
        character(len=*), intent(in) :: output_filename

        type(JSON_Value_Pointer) :: output_data
        type(JSON_Value_Pointer), allocatable :: json_particles(:)
        integer :: i_particle
        character(len=1024) :: string_i

        call json_create_object(output_data%ptr, "")
        allocate(json_particles(particles%number%get()))
        do i_particle = 1, particles%number%get()
            json_particles(i_particle)%ptr => null()
            write(string_i, *) i_particle
            call json_create_object(json_particles(i_particle)%ptr, &
                                    "Particle "//trim(adjustl(string_i)))
            call json_add(output_data%ptr, json_particles(i_particle)%ptr)
        end do

        do i_particle = 1, particles%number%get()
            call json_add(json_particles(i_particle)%ptr, &
                "diameter", particles%diameter%get())
            call json_add(json_particles(i_particle)%ptr, &
                "moment norm", particles%moment_norm%get())
            call json_add(json_particles(i_particle)%ptr, &
                "position", particles%positions%get(i_particle))
            call json_add(json_particles(i_particle)%ptr, &
                "orientation", particles%orientations%get(i_particle))
        end do

        call json_print(output_data%ptr, output_filename)
        deallocate(json_particles)
        call json_destroy(output_data%ptr)
    end subroutine json_write_particles

end module procedures_particles_exchange_write

program test_particles_exchange

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_geometry, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use procedures_random, only: random_integer, random_orientation
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box
use types_particle, only: Concrete_Particle
use types_particles, only: Particles_Wrapper
use class_particles_factory, only: Abstract_Particles_Factory, &
    Dipolar_Particles_Factory, Apolar_Particles_Factory
use class_particles_exchange, only: Particles_Exchange_Facade
use procedures_particles_exchange_write, only: json_write_particles

implicit none

    type(Particles_Exchange_Facade) :: particles_exchange
    class(Abstract_Particles_Factory), allocatable :: particles_factory
    type(Particles_Wrapper) :: particles
    type(Concrete_Particle) :: particle
    class(Abstract_Periodic_Box), allocatable :: periodic_box
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found

    real(DP), allocatable :: box_size(:)
    logical :: dipolar
    real(DP), dimension(num_dimensions) :: rand_3d
    integer :: i_particle

    call json_initialize()
    data_filename = "particles_exchange.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    data_field = "Periodic Box.size"
    call input_data%get(data_field, box_size, data_found)
    call test_data_found(data_field, data_found)
    allocate(XYZ_Periodic_Box :: periodic_box)
    call periodic_box%set(box_size)

    data_field = "Particles.dipolar"
    call input_data%get(data_field, dipolar, data_found)
    call test_data_found(data_field, data_found)
    if (dipolar) then
        allocate(Dipolar_Particles_Factory :: particles_factory)
    else
        allocate(Apolar_Particles_Factory :: particles_factory)
    end if
    call particles_factory%create(particles, periodic_box, input_data, "Particles")

    call particles_exchange%construct(particles)
    call json_write_particles(particles, "initial.json")

    call random_number(rand_3d)
    particle%position = periodic_box%get_size() * rand_3d
    particle%orientation = random_orientation()
    call particles_exchange%add(particle)
    call json_write_particles(particles, "added.json")

    i_particle = random_integer(particles%number%get())
    write(output_unit, *) "remove index =", i_particle
    call particles_exchange%remove(i_particle)
    call json_write_particles(particles, "removed.json")

    call particles_exchange%destroy()
    call particles_factory%destroy(particles)
    deallocate(particles_factory)
    deallocate(periodic_box)
    deallocate(data_field)
    call input_data%destroy()

end program test_particles_exchange
