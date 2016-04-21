module procedures_des_reci_factory

use procedures_des_reci_weight_factory, only: des_reci_weight_create => create, &
    des_reci_weight_destroy => destroy
use procedures_des_reci_structure_factory, only: des_reci_structure_create => create, &
    des_reci_structure_destroy => destroy
use procedures_des_reci_visitor_factory, only: des_reci_visitor_create => create, &
    des_reci_visitor_destroy => destroy

implicit none

private
public :: create, destroy

interface create
    module procedure :: des_reci_visitor_create
    module procedure :: des_reci_structure_create
    module procedure :: des_reci_weight_create
end interface create

interface destroy
    module procedure :: des_reci_weight_destroy
    module procedure :: des_reci_structure_destroy
    module procedure :: des_reci_visitor_destroy
end interface destroy

end module procedures_des_reci_factory
