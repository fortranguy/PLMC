program test_canonical

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use types_box, only: Box_Wrapper
use procedures_box_factory, only: box_factory_create, box_factory_destroy
use types_particles, only: Mixture_Wrapper
use procedures_particles_factory, only: mixture_factory_create, mixture_factory_destroy
use types_changes, only: Changes_Wrapper
use procedures_changes_factory, only: changes_factory_create, changes_factory_destroy
use types_short_potential, only: Short_Potential_Wrapper
use procedures_short_potential_factory, only: short_potential_factory_create, &
    short_potential_factory_destroy
use class_one_particle_move, only: Metropolis_One_Particle_Move

implicit none

    type(Box_Wrapper) :: box
    type(Mixture_Wrapper) :: mixture
    type(Changes_Wrapper) :: changes
    type(Short_Potential_Wrapper) :: short_potential
    type(Metropolis_One_Particle_Move) :: one_particle_move

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, prefix
    logical :: data_found
    real(DP), allocatable :: box_size(:)
    logical :: dipolar

    call json_initialize()
    data_filename = "canonical.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)
    prefix = "Test Canonical"

    call box_factory_create(box, input_data, prefix)

    call box_factory_destroy(box)
    deallocate(prefix)
    deallocate(data_filename)
    call input_data%destroy()

end program test_canonical
