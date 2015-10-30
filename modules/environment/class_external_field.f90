module class_external_field

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use class_field_expression, only: Abstract_Field_Expression
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain

implicit none

private

    type, abstract, public :: Abstract_External_Field
    private
        class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
        class(Abstract_Field_Expression), allocatable :: field_expression
    contains
        procedure :: construct => Abstract_External_Field_construct
        procedure :: destroy => Abstract_External_Field_destroy
        procedure :: get => Abstract_External_Field_get
    end type Abstract_External_Field

    type, extends(Abstract_External_Field), public :: Concrete_External_Field

    end type Concrete_External_Field

    type, extends(Abstract_External_Field), public :: Null_External_Field
    contains
        procedure :: construct => Null_External_Field_construct
        procedure :: destroy => Null_External_Field_destroy
        procedure :: get => Null_External_Field_get
    end type Null_External_Field

contains

!implementation Abstract_External_Field

    subroutine Abstract_External_Field_construct(this, parallelepiped_domain, field_expression)
        class(Abstract_External_Field), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), intent(in) :: parallelepiped_domain
        class(Abstract_Field_Expression), intent(in) :: field_expression

        allocate(this%parallelepiped_domain, source = parallelepiped_domain)
        allocate(this%field_expression, source = field_expression)
    end subroutine Abstract_External_Field_construct

    subroutine Abstract_External_Field_destroy(this)
        class(Abstract_External_Field), intent(inout) :: this

        if (allocated(this%field_expression)) deallocate(this%field_expression)
        if (allocated(this%parallelepiped_domain)) deallocate(this%parallelepiped_domain)
    end subroutine Abstract_External_Field_destroy

    pure function Abstract_External_Field_get(this, position) result(external_field)
        class(Abstract_External_Field), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: external_field(num_dimensions)

        if (this%parallelepiped_domain%is_inside(position)) then
            external_field = this%field_expression%get(position)
        else
            external_field = 0._DP
        end if
    end function Abstract_External_Field_get

!end implementation Abstract_External_Field

!implementation Null_External_Field

    subroutine Null_External_Field_construct(this, parallelepiped_domain, field_expression)
        class(Null_External_Field), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), intent(in) :: parallelepiped_domain
        class(Abstract_Field_Expression), intent(in) :: field_expression
    end subroutine Null_External_Field_construct

    subroutine Null_External_Field_destroy(this)
        class(Null_External_Field), intent(inout) :: this
    end subroutine Null_External_Field_destroy

    pure function Null_External_Field_get(this, position) result(external_field)
        class(Null_External_Field), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: external_field(num_dimensions)
        external_field = 0._DP
    end function Null_External_Field_get

!end implementation Null_External_Field

end module class_external_field