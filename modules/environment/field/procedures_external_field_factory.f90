module procedures_external_field_factory

use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_field_expression, only: Abstract_Field_Expression
use classes_external_field, only: Abstract_External_Field, Concrete_External_Field, &
    Null_External_Field

implicit none

private
public :: external_field_create, external_field_destroy

contains

    subroutine external_field_create(external_field, parallelepiped_domain, field_expression, &
        field_applied)
        class(Abstract_External_Field), allocatable, intent(out) :: external_field
        class(Abstract_Parallelepiped_Domain), intent(in) :: parallelepiped_domain
        class(Abstract_Field_Expression), intent(in) :: field_expression
        logical, intent(in) :: field_applied

        if (field_applied) then
            allocate(Concrete_External_Field :: external_field)
        else
            allocate(Null_External_Field :: external_field)
        end if
        call external_field%construct(parallelepiped_domain, field_expression)
    end subroutine external_field_create

    subroutine external_field_destroy(external_field)
        class(Abstract_External_Field), allocatable, intent(inout) :: external_field

        call external_field%destroy()
        if (allocated(external_field)) deallocate(external_field)
    end subroutine external_field_destroy

end module procedures_external_field_factory
