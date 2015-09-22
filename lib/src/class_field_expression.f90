module class_field_expression

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use types_field_parameters, only: Abstract_Field_Parameters, Constant_Field_Parameters, &
    Null_Field_Parameters
use procedures_errors, only: error_exit

implicit none

private

    type, abstract, public :: Abstract_Field_Expression
    contains
        procedure(Abstract_Field_Expression_set), deferred :: set
        procedure(Abstract_Field_Expression_get), deferred :: get
    end type Abstract_Field_Expression

    abstract interface

        subroutine Abstract_Field_Expression_set(this, parameters)
        import :: Abstract_Field_Parameters, Abstract_Field_Expression
            class(Abstract_Field_Expression), intent(inout) :: this
            class(Abstract_Field_Parameters), intent(in) :: parameters
        end subroutine Abstract_Field_Expression_set

        pure function Abstract_Field_Expression_get(this, position) result(expression)
        import :: DP, num_dimensions, Abstract_Field_Expression
            class(Abstract_Field_Expression), intent(in) :: this
            real(DP), intent(in) :: position(:)
            real(DP) :: expression(num_dimensions)
        end function Abstract_Field_Expression_get

    end interface

    type, extends(Abstract_Field_Expression), public :: Null_Field_Expression
    private
        type(Constant_Field_Parameters) :: parameters
    contains
        procedure :: set => Null_Field_Expression_set
        procedure :: get => Null_Field_Expression_get
    end type Null_Field_Expression

    type, extends(Abstract_Field_Expression), public :: Constant_Field_Expression
    private
        type(Constant_Field_Parameters) :: parameters
    contains
        procedure :: set => Constant_Field_Expression_set
        procedure :: get => Constant_Field_Expression_get
    end type Constant_Field_Expression

contains

!implementation Constant_Field_Expression_get

    subroutine Constant_Field_Expression_set(this, parameters)
        class(Constant_Field_Expression), intent(inout) :: this
        class(Abstract_Field_Parameters), intent(in) :: parameters

        select type (parameters)
            type is (Constant_Field_Parameters)
                this%parameters%vector = parameters%vector
            class default
                call error_exit("Constant_Field_Expression: no parameters were given.")
        end select
    end subroutine Constant_Field_Expression_set

    pure function Constant_Field_Expression_get(this, position) result(expression)
        class(Constant_Field_Expression), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: expression(num_dimensions)

        expression = this%parameters%vector
    end function Constant_Field_Expression_get

!end implementation Constant_Field_Expression_get

!implementation Null_Field_Expression_get

    subroutine Null_Field_Expression_set(this, parameters)
        class(Null_Field_Expression), intent(inout) :: this
        class(Abstract_Field_Parameters), intent(in) :: parameters
    end subroutine Null_Field_Expression_set

    pure function Null_Field_Expression_get(this, position) result(expression)
        class(Null_Field_Expression), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: expression(num_dimensions)
        expression = 0._DP
    end function Null_Field_Expression_get

!end implementation Null_Field_Expression_get

end module class_field_expression
