program test_canonical

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use types_box, only: Box_Wrapper
use procedures_box_factory, only: box_factory_create, box_factory_destroy
use types_particles, only: Mixture_Wrapper
use procedures_particles_factory, only: particles_factory_create, particles_factory_destroy
use types_changes, only: Changes_Wrapper
use procedures_changes_factory, only: changes_factory_create, changes_factory_destroy
use types_short_potential, only: Mixture_Short_Potentials_Wrapper
use procedures_short_potential_factory, only: short_potential_factory_create, &
    short_potential_factory_destroy
use class_one_particle_move, only: Abstract_One_Particle_Move

implicit none

    type(Box_Wrapper) :: box
    type(Mixture_Wrapper) :: mixture
    type(Changes_Wrapper) :: changes_1, changes_2
    type(Mixture_Short_Potentials_Wrapper) :: mixture_short_potentials
    class(Abstract_One_Particle_Move), allocatable :: one_particle_move

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()
    data_filename = "canonical.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call box_factory_create(box, input_data, "Box.")
    call particles_factory_create(mixture%components(1), input_data, "Mixture.Component 1.", &
        box%periodic_box)
    call particles_factory_create(mixture%components(2), input_data, "Mixture.Component 2.", &
        box%periodic_box)
    call particles_factory_create(mixture%inter_diameter, mixture%components(1)%diameter, &
        mixture%components(2)%diameter, input_data, "Mixture.Inter 12.")
    call changes_factory_create(changes_1, input_data, "Changes.Component 1.", &
        mixture%components(1))
    call changes_factory_create(changes_2, input_data, "Changes.Component 2.", &
        mixture%components(2))
    call short_potential_factory_create(mixture_short_potentials%intras(1), input_data, &
        "Short Potentials.Component 1.", box%periodic_box, mixture%components(1))
    call short_potential_factory_create(mixture_short_potentials%intras(2), input_data, &
        "Short Potentials.Component 2.", box%periodic_box, mixture%components(2))
    call short_potential_factory_create(mixture_short_potentials%inter_micro, input_data, &
        "Short Potentials.Inter 12.", mixture%inter_diameter)
    call short_potential_factory_create(mixture_short_potentials%inter_macros(1), &
        mixture_short_potentials%inter_micro, input_data, "Short Potentials.Inter 12.", &
        box%periodic_box, mixture%components(1)%positions)
    call short_potential_factory_create(mixture_short_potentials%inter_macros(2), &
        mixture_short_potentials%inter_micro, input_data, "Short Potentials.Inter 12.", &
        box%periodic_box, mixture%components(2)%positions)
    call input_data%destroy()

    call short_potential_factory_destroy(mixture_short_potentials%inter_macros(2))
    call short_potential_factory_destroy(mixture_short_potentials%inter_macros(1))
    call short_potential_factory_destroy(mixture_short_potentials%inter_micro)
    call short_potential_factory_destroy(mixture_short_potentials%intras(2))
    call short_potential_factory_destroy(mixture_short_potentials%intras(1))
    call changes_factory_destroy(changes_1)
    call changes_factory_destroy(changes_2)
    call particles_factory_destroy(mixture%inter_diameter)
    call particles_factory_destroy(mixture%components(2))
    call particles_factory_destroy(mixture%components(1))
    call box_factory_destroy(box)

end program test_canonical
