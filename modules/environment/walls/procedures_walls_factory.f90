module procedures_walls_factory

use procedures_floor_penetration_factory, only: floor_penetration_create => create, &
    floor_penetration_destroy => destroy
use procedures_walls_potential_factory, only: walls_potential_create => create, &
    walls_potential_destroy => destroy
use procedures_walls_potential_visitor_factory, only: walls_potential_visitor_create => create, &
    walls_potential_visitor_destroy => destroy

implicit none

private
public :: create, destroy

interface create
    module procedure :: floor_penetration_create
    module procedure :: walls_potential_create
    module procedure :: walls_potential_visitor_create
end interface create

interface destroy
    module procedure :: walls_potential_visitor_destroy
    module procedure :: walls_potential_destroy
    module procedure :: floor_penetration_destroy
end interface destroy

end module procedures_walls_factory
