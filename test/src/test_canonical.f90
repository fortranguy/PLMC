program test_canonical

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box
use module_particles, only: Particles_Wrapper
use class_particles_factory, only: Concrete_Particles_Factory
use class_one_particle_move, only: Metropolis_One_Particle_Move

implicit none

    type(Concrete_Particles_Factory) :: particles_factory
    type(Particles_Wrapper) :: particles
    class(Abstract_Periodic_Box), allocatable :: periodic_box
    type(Metropolis_One_Particle_Move) :: one_particle_move
    type(json_file) :: input_data

    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found
    real(DP), allocatable :: box_size(:)
    logical :: dipolar

    call json_initialize()
    data_filename = "canonical.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call input_data%destroy()

end program test_canonical
