module procedures_box_factory

use procedures_periodic_box_factory, only: periodic_box_create => create, &
    periodic_box_destroy => destroy
use procedures_box_size_memento_factory, only: box_size_memento_create => create, &
    box_size_memento_destroy => destroy
use procedures_parallelepiped_domain_factory, only: parallelepiped_domain_create_from_json => &
    create_from_json, parallelepiped_domain_create_from_box => create_from_box, &
    parallelepiped_domain_create_from_walls => create_from_walls, parallelepiped_domain_destroy => &
    destroy
use procedures_reciprocal_lattice_factory, only: reciprocal_lattice_create => create, &
    reciprocal_lattice_destroy => destroy
use procedures_box_size_checker_factory, only: box_size_checker_create => create, &
    box_size_checker_destroy => destroy

implicit none

private
public :: create, destroy

interface create
    module procedure :: periodic_box_create
    module procedure :: box_size_memento_create
    module procedure :: parallelepiped_domain_create_from_json
    module procedure :: parallelepiped_domain_create_from_box
    module procedure :: parallelepiped_domain_create_from_walls
    module procedure :: reciprocal_lattice_create
    module procedure :: box_size_checker_create
end interface create

interface destroy
    module procedure :: box_size_checker_destroy
    module procedure :: reciprocal_lattice_destroy
    module procedure :: parallelepiped_domain_destroy
    module procedure :: box_size_memento_destroy
    module procedure :: periodic_box_destroy
end interface destroy

end module procedures_box_factory
