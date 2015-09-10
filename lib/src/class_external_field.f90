module class_external_field

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_field_expression, only: Abstract_Field_Expression
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain

implicit none

private

    type, public :: External_Field_Facade
        class(Abstract_Parallelepiped_Domain), pointer :: parallelepiped_domain
        class(Abstract_Field_Expression), pointer :: field_expression
    contains
        procedure :: construct => External_Field_Facade_construct
        procedure :: destroy => External_Field_Facade_destroy
        procedure :: get => External_Field_Facade_get
    end type External_Field_Facade
    
contains

!implementation External_Field_Facade

    subroutine External_Field_Facade_construct(this, parallelepiped_domain, field_expression)
        class(External_Field_Facade), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: parallelepiped_domain
        class(Abstract_Field_Expression), target, intent(in) :: field_expression
        
        this%parallelepiped_domain => parallelepiped_domain
        this%field_expression => field_expression
    end subroutine External_Field_Facade_construct

    subroutine External_Field_Facade_destroy(this)
        class(External_Field_Facade), intent(inout) :: this

        this%field_expression => null()
        this%parallelepiped_domain => null()
    end subroutine External_Field_Facade_destroy

    pure function External_Field_Facade_get(this, position) result(external_field)
        class(External_Field_Facade), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: external_field(num_dimensions)

        if (this%parallelepiped_domain%is_inside(position)) then
            external_field = this%field_expression%get(position)
        else
            external_field = 0._DP
        end if
    end function External_Field_Facade_get

!end implementation External_Field_Facade
    
end module class_external_field
