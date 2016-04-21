module procedures_field_factory

use procedures_field_expression_factory, only: field_expression_create, field_expression_destroy
use procedures_external_field_factory, only: external_field_create, external_field_destroy

implicit none

private
public :: field_create, field_destroy

interface field_create
    module procedure :: field_expression_create
    module procedure :: external_field_create
end interface field_create

interface field_destroy
    module procedure :: external_field_destroy
    module procedure :: field_expression_destroy
end interface field_destroy

end module procedures_field_factory
