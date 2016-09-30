module procedures_external_fields_factory

use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_field_expression, only: Abstract_Field_Expression
use classes_external_field, only: Abstract_External_Field, Concrete_External_Field, &
    Null_External_Field

implicit none

private
public :: create, destroy

contains

    subroutine create(external_fields, parallelepiped_domains, field_expression, field_applied)
        class(Abstract_External_Field), allocatable, intent(out) :: external_fields(:)
        class(Abstract_Parallelepiped_Domain), intent(in) :: parallelepiped_domains(:)
        class(Abstract_Field_Expression), intent(in) :: field_expression
        logical, intent(in) :: field_applied

        integer :: i_box

        if (field_applied) then
            allocate(Concrete_External_Field :: external_fields(size(parallelepiped_domains)))
        else
            allocate(Null_External_Field :: external_fields(size(parallelepiped_domains)))
        end if
        do i_box = 1, size(external_fields)
            call external_fields(i_box)%construct(parallelepiped_domains(i_box), field_expression)
        end do
    end subroutine create

    subroutine destroy(external_fields)
        class(Abstract_External_Field), allocatable, intent(inout) :: external_fields(:)

        integer :: i_box

        if (allocated(external_fields)) then
            do i_box = size(external_fields), 1, -1
                call external_fields(i_box)%destroy()
            end do
            deallocate(external_fields)
        end if
    end subroutine destroy

end module procedures_external_fields_factory
