module procedures_walls_factory

use procedures_floor_penetration_factory, only: floor_penetration_create => create, &
    floor_penetration_destroy => destroy
use procedures_visitable_walls_factory, only: visitable_walls_create => create, &
    visitable_walls_destroy => destroy
use procedures_walls_visitor_factory, only: walls_visitor_create => create, &
    walls_visitor_destroy => destroy

implicit none

private
public :: create, destroy

interface create
    module procedure :: floor_penetration_create
    module procedure :: visitable_walls_create
    module procedure :: walls_visitor_create
end interface create

interface destroy
    module procedure :: walls_visitor_destroy
    module procedure :: visitable_walls_destroy
    module procedure :: floor_penetration_destroy
end interface destroy

end module procedures_walls_factory
