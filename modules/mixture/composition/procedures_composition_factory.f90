module procedures_composition_factory

use procedures_num_particles_factory, only: num_particles_create_line => create_line, &
    num_particles_create_element => create_element, num_particles_detroy_line => destroy_line, &
    num_particles_detroy_element => destroy_element
use procedures_component_chemical_potential_factory, only: &
    component_chemical_potential_create => create, component_chemical_potential_destroy => destroy
use procedures_average_num_particles_factory, only: average_num_particles_create => create, &
    average_num_particles_destroy => destroy

implicit none

private
public :: create, destroy

interface create
    module procedure :: num_particles_create_line, num_particles_create_element
    module procedure :: component_chemical_potential_create
    module procedure :: average_num_particles_create
end interface create

interface destroy
    module procedure :: average_num_particles_destroy
    module procedure :: component_chemical_potential_destroy
    module procedure :: num_particles_detroy_line, num_particles_detroy_element
end interface destroy

end module procedures_composition_factory
