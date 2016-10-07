module procedures_composition_factory

use procedures_num_particles_factory, only: num_particles_create => create, &
    num_particles_detroy => destroy
use procedures_component_chemical_potential_factory, only: &
    component_chemical_potential_create => create, component_chemical_potential_destroy => destroy
use procedures_component_average_number_factory, only: component_average_number_create => create, &
    component_average_number_destroy => destroy

implicit none

private
public :: create, destroy

interface create
    module procedure :: num_particles_create
    module procedure :: component_chemical_potential_create
    module procedure :: component_average_number_create
end interface create

interface destroy
    module procedure :: component_average_number_destroy
    module procedure :: component_chemical_potential_destroy
    module procedure :: num_particles_detroy
end interface destroy

end module procedures_composition_factory
