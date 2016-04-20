module procedures_des_reci_factory

use procedures_des_reci_weight_factory, only: des_reci_weight_create, des_reci_weight_destroy
use procedures_des_reci_structure_factory, only: des_reci_structure_create, &
    des_reci_structure_destroy
use procedures_des_reci_visitor_factory, only: des_reci_visitor_create, des_reci_visitor_destroy

implicit none

private
public :: des_reci_create, des_reci_destroy

interface des_reci_create
    module procedure :: des_reci_visitor_create
    module procedure :: des_reci_structure_create
    module procedure :: des_reci_weight_create
end interface des_reci_create

interface des_reci_destroy
    module procedure :: des_reci_weight_destroy
    module procedure :: des_reci_structure_destroy
    module procedure :: des_reci_visitor_destroy
end interface des_reci_destroy

end module procedures_des_reci_factory
