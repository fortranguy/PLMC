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
use module_particles, only: Concrete_Particles

implicit none

private
public json_write_particles

contains

    subroutine json_write_particles(particles, output_filename)
        type(Concrete_Particles), intent(in) :: particles
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
use procedures_orientation, only: random_orientation
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box
use class_number, only: Abstract_Number, Concrete_Number
use class_diameter, only: Abstract_Diameter, Concrete_Diameter
use class_moment_norm, only: Abstract_Moment_Norm, Concrete_Moment_Norm
use class_positions, only: Abstract_Positions, Concrete_Positions
use class_orientations, only: Abstract_Orientations, Concrete_Orientations
use module_particles, only: Concrete_Particle, Concrete_Particles, &
    Concrete_Particles_construct, Concrete_Particles_destroy
use class_particles_exchange, only: Particles_Exchange_Facade
use procedures_particles_exchange_write, only: json_write_particles

implicit none

    type(Particles_Exchange_Facade) :: particles_exchange
    type(Concrete_Particles) :: particles
    type(Concrete_Particle) :: particle
    class(Abstract_Periodic_Box), allocatable :: periodic_box
    class(Abstract_Number), allocatable :: number
    class(Abstract_Diameter), allocatable :: diameter
    class(Abstract_Moment_Norm), allocatable :: moment_norm
    class(Abstract_Positions), allocatable :: positions
    class(Abstract_Orientations), allocatable :: orientations

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found

    real(DP), allocatable :: box_size(:)
    real(DP), dimension(num_dimensions) :: rand_3d
    real(DP) :: diameter_value, diameter_min_factor, moment_norm_value, rand
    integer :: num_particles, i_particle

    call json_initialize()
    data_filename = "particles_exchange.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    data_field = "Periodic Box.size"
    call input_data%get(data_field, box_size, data_found)
    call test_data_found(data_field, data_found)
    allocate(XYZ_Periodic_Box :: periodic_box)
    call periodic_box%set_size(box_size)

    data_field = "Particles.number"
    call input_data%get(data_field, num_particles, data_found)
    call test_data_found(data_field, data_found)
    allocate(Concrete_Number :: number)
    call number%set(num_particles)

    allocate(Concrete_Diameter :: diameter)
    data_field = "Particles.diameter"
    call input_data%get(data_field, diameter_value, data_found)
    call test_data_found(data_field, data_found)
    data_field = "Particles.diameter minimum factor"
    call input_data%get(data_field, diameter_min_factor, data_found)
    call test_data_found(data_field, data_found)
    call diameter%set(diameter_value, diameter_min_factor)

    data_field = "Particles.moment norm"
    call input_data%get(data_field, moment_norm_value, data_found)
    call test_data_found(data_field, data_found)
    allocate(Concrete_Moment_Norm :: moment_norm)
    call moment_norm%set(moment_norm_value)

    allocate(Concrete_Positions :: positions)
    call positions%construct(periodic_box, number)
    do i_particle = 1, positions%get_num()
        call random_number(rand_3d)
        call positions%set(i_particle, periodic_box%get_size() * rand_3d)
    end do

    allocate(Concrete_Orientations :: orientations)
    call orientations%construct(number)
    do i_particle = 1, orientations%get_num()
        call orientations%set(i_particle, random_orientation())
    end do

    call Concrete_Particles_construct(particles, number, diameter, moment_norm, &
                                     positions, orientations)
    call particles_exchange%construct(particles)
    call json_write_particles(particles, "initial.json")

    particle%diameter = diameter%get()
    particle%min_diameter = diameter%get_min()
    particle%moment_norm = moment_norm%get()
    call random_number(rand_3d)
    particle%position = periodic_box%get_size() * rand_3d
    particle%orientation = random_orientation()
    call particles_exchange%add(particle)
    call json_write_particles(particles, "added.json")

    call random_number(rand)
    i_particle = int(real(number%get(), DP) * rand) + 1
    write(output_unit, *) "remove index =", i_particle
    call particles_exchange%remove(i_particle)
    call json_write_particles(particles, "removed.json")

    call particles_exchange%destroy()
    call Concrete_Particles_destroy(particles)
    call orientations%destroy()
    deallocate(orientations)
    call positions%destroy()
    deallocate(positions)
    deallocate(moment_norm)
    deallocate(diameter)
    deallocate(number)
    deallocate(periodic_box)
    deallocate(data_field)

    call input_data%destroy()

end program test_particles_exchange
