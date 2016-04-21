module procedures_walls_factory

use procedures_floor_penetration_factory, only: floor_penetration_create, floor_penetration_destroy
use procedures_walls_potential_factory, only: walls_potential_create, walls_potential_destroy
use procedures_walls_potential_visitor_factory, only: walls_potential_visitor_create, &
    walls_potential_visitor_destroy

implicit none

private
public :: walls_create, walls_destroy

interface walls_create
    module procedure :: floor_penetration_create
    module procedure :: walls_potential_create
    module procedure :: walls_potential_visitor_create
end interface walls_create

interface walls_destroy
    module procedure :: walls_potential_visitor_destroy
    module procedure :: walls_potential_destroy
    module procedure :: floor_penetration_destroy
end interface walls_destroy

end module procedures_walls_factory
