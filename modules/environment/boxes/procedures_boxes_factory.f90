module procedures_boxes_factory

use procedures_periodic_boxes_factory, only: periodic_boxes_create => create, &
    periodic_boxes_destroy => destroy
use procedures_box_size_memento_factory, only: box_size_memento_create => create, &
    box_size_memento_destroy => destroy
use procedures_parallelepiped_domains_factory, only: parallelepiped_domains_create_from_json => &
    create_from_json, parallelepiped_domains_create_from_boxes => create_from_boxes, &
    parallelepiped_domains_create_from_walls => create_from_walls, parallelepiped_domains_destroy => &
    destroy
use procedures_reciprocal_lattices_factory, only: reciprocal_lattices_create => create, &
    reciprocal_lattices_destroy => destroy
use procedures_box_size_checkers_factory, only: box_size_checkers_create => create, &
    box_size_checkers_destroy => destroy

implicit none

private
public :: create, destroy

interface create
    module procedure :: periodic_boxes_create
    module procedure :: box_size_memento_create
    module procedure :: parallelepiped_domains_create_from_json
    module procedure :: parallelepiped_domains_create_from_boxes
    module procedure :: parallelepiped_domains_create_from_walls
    module procedure :: reciprocal_lattices_create
    module procedure :: box_size_checkers_create
end interface create

interface destroy
    module procedure :: box_size_checkers_destroy
    module procedure :: reciprocal_lattices_destroy
    module procedure :: parallelepiped_domains_destroy
    module procedure :: box_size_memento_destroy
    module procedure :: periodic_boxes_destroy
end interface destroy

end module procedures_boxes_factory
