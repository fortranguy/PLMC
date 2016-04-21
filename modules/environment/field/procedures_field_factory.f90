module procedures_field_factory

use procedures_field_expression_factory, only: field_expression_create => create, &
    field_expression_destroy => destroy
use procedures_external_field_factory, only: external_field_create => create, &
    external_field_destroy => destroy

implicit none

private
public :: create, destroy

interface create
    module procedure :: field_expression_create
    module procedure :: external_field_create
end interface create

interface destroy
    module procedure :: external_field_destroy
    module procedure :: field_expression_destroy
end interface destroy

end module procedures_field_factory
