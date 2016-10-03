module procedures_dlc_factory

use procedures_dlc_weight_factory, only: dlc_weight_create => create, dlc_weight_destroy => destroy
use procedures_dlc_structures_factory, only: dlc_structures_create => create, &
    dlc_structures_destroy => destroy
use procedures_dlc_visitors_factory, only: dlc_visitors_create => create, &
    dlc_visitors_destroy => destroy

implicit none

private
public :: create, destroy

interface create
    module procedure :: dlc_visitors_create
    module procedure :: dlc_structures_create
    module procedure :: dlc_weight_create
end interface create

interface destroy
    module procedure :: dlc_weight_destroy
    module procedure :: dlc_structures_destroy
    module procedure :: dlc_visitors_destroy
end interface destroy

end module procedures_dlc_factory
