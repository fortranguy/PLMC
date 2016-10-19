module classes_external_field

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use classes_field_expression, only: Abstract_Field_Expression
use procedures_field_expression_factory, only: field_expression_destroy => destroy
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain

implicit none

private

    type, abstract, public :: Abstract_External_Field
    private
        class(Abstract_Parallelepiped_Domain), pointer :: field_domain => null()
        class(Abstract_Field_Expression), allocatable :: field_expression
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: get => Abstract_get
    end type Abstract_External_Field

    type, extends(Abstract_External_Field), public :: Concrete_External_Field

    end type Concrete_External_Field

    type, extends(Abstract_External_Field), public :: Null_External_Field
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: get => Null_get
    end type Null_External_Field

contains

!implementation Abstract_External_Field

    subroutine Abstract_construct(this, field_domain, field_expression)
        class(Abstract_External_Field), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: field_domain
        class(Abstract_Field_Expression), intent(in) :: field_expression

        this%field_domain => field_domain
        allocate(this%field_expression, source = field_expression)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_External_Field), intent(inout) :: this

        call field_expression_destroy(this%field_expression)
        this%field_domain => null()
    end subroutine Abstract_destroy

    pure function Abstract_get(this, position) result(external_field)
        class(Abstract_External_Field), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: external_field(num_dimensions)

        if (this%field_domain%is_inside(position)) then
            external_field = this%field_expression%get(position)
        else
            external_field = 0._DP
        end if
    end function Abstract_get

!end implementation Abstract_External_Field

!implementation Null_External_Field

    subroutine Null_construct(this, field_domain, field_expression)
        class(Null_External_Field), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: field_domain
        class(Abstract_Field_Expression), intent(in) :: field_expression
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_External_Field), intent(inout) :: this
    end subroutine Null_destroy

    pure function Null_get(this, position) result(external_field)
        class(Null_External_Field), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: external_field(num_dimensions)
        external_field = 0._DP
    end function Null_get

!end implementation Null_External_Field

end module classes_external_field
