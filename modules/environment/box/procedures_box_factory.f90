module procedures_box_factory

use procedures_periodic_box_factory, only: periodic_box_create, periodic_box_destroy
use procedures_parallelepiped_domain_factory, only: parallelepiped_domain_create, &
    parallelepiped_domain_destroy
use procedures_reciprocal_lattice_factory, only: reciprocal_lattice_create, &
    reciprocal_lattice_destroy
use procedures_box_size_checker_factory, only: box_size_checker_create, box_size_checker_destroy

implicit none

private
public :: box_create, box_destroy

interface box_create
    module procedure :: periodic_box_create
    module procedure :: parallelepiped_domain_create
    module procedure :: reciprocal_lattice_create
    module procedure :: box_size_checker_create
end interface box_create

interface box_destroy
    module procedure :: box_size_checker_destroy
    module procedure :: reciprocal_lattice_destroy
    module procedure :: parallelepiped_domain_destroy
    module procedure :: periodic_box_destroy
end interface box_destroy

end module procedures_box_factory
