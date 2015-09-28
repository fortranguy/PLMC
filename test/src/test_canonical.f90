module types_canonical

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

    type, public :: Concrete_Particles_Observables
        integer :: num_move_hits = 0
        integer :: num_move_rejection = 0
        real(DP) :: energy_step = 0._DP
    end type Concrete_Particles_Observables

end module types_canonical

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
use class_one_particle_move, only: Abstract_One_Particle_Move, &
    Two_Candidates_One_Particle_Move
use types_canonical, only: Concrete_Particles_Observables

implicit none

    type(Box_Wrapper) :: box
    type(Mixture_Wrapper) :: mixture
    type(Changes_Wrapper) :: changes_1, changes_2
    type(Mixture_Short_Potentials_Wrapper) :: mixture_short_potentials
    class(Abstract_One_Particle_Move), allocatable :: one_particle_move
    type(Concrete_Particles_Observables) :: observable_1, observable_2
    logical :: move_success

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
    allocate(Two_Candidates_One_Particle_Move :: one_particle_move)
    call one_particle_move%construct(box%temperature)
    call one_particle_move%set_first_candidate(mixture%components(1)%positions, &
        changes_1%moved_positions, mixture_short_potentials%intras(1)%cells, &
        mixture_short_potentials%inter_macros(1)%cells)
    call one_particle_move%set_second_candidate(mixture%components(2)%positions, &
        changes_2%moved_positions, mixture_short_potentials%intras(2)%cells, &
        mixture_short_potentials%inter_macros(2)%cells)

    !call one_particle_move%try(move_success, )

    deallocate(one_particle_move)
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
