module procedures_dlc_factory

use procedures_dlc_weight_factory, only: dlc_weight_create, dlc_weight_destroy
use procedures_dlc_structures_factory, only: dlc_structures_create, dlc_structures_destroy
use procedures_dlc_visitor_factory, only: dlc_visitor_create, dlc_visitor_destroy

implicit none

private
public :: dlc_create, dlc_destroy

interface dlc_create
    module procedure :: dlc_visitor_create
    module procedure :: dlc_structures_create
    module procedure :: dlc_weight_create
end interface dlc_create

interface dlc_destroy
    module procedure :: dlc_weight_destroy
    module procedure :: dlc_structures_destroy
    module procedure :: dlc_visitor_destroy
end interface dlc_destroy

end module procedures_dlc_factory
