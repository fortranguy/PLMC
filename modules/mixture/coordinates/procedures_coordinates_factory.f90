module procedures_coordinates_factory

use procedures_component_coordinates_factory, only: &
    component_coordinates_create_positions => create_positions, &
    component_coordinates_create_orientations => create_orientations, &
    component_coordinates_destroy => destroy
use procedures_component_dipole_moments_factory, only: component_dipole_moments_create => create,&
    component_dipole_moments_destroy => destroy

implicit none

private
public :: create, destroy

interface create
    module procedure :: component_coordinates_create_positions
    module procedure :: component_coordinates_create_orientations
    module procedure :: component_dipole_moments_create
end interface create

interface destroy
    module procedure :: component_dipole_moments_destroy
    module procedure :: component_coordinates_destroy
end interface destroy

end module procedures_coordinates_factory
