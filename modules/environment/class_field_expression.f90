module class_field_expression

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_3d_array

implicit none

private

    type, abstract, public :: Abstract_Field_Expression
    contains
        procedure(Abstract_get), deferred :: get
    end type Abstract_Field_Expression

    abstract interface

        pure function Abstract_get(this, position) result(expression)
        import :: DP, num_dimensions, Abstract_Field_Expression
            class(Abstract_Field_Expression), intent(in) :: this
            real(DP), intent(in) :: position(:)
            real(DP) :: expression(num_dimensions)
        end function Abstract_get

    end interface

    type, extends(Abstract_Field_Expression), public :: Null_Field_Expression
    contains
        procedure :: set => Null_set
        procedure :: get => Null_get
    end type Null_Field_Expression

    type, extends(Abstract_Field_Expression), public :: Constant_Field_Expression
    private
        real(DP) :: vector(num_dimensions)
    contains
        procedure :: set => Constant_set
        procedure :: get => Constant_get
    end type Constant_Field_Expression

contains

!implementation Null_Field_Expression

    subroutine Null_set(this)
        class(Null_Field_Expression), intent(inout) :: this
    end subroutine Null_set

    pure function Null_get(this, position) result(expression)
        class(Null_Field_Expression), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: expression(num_dimensions)
        expression = 0._DP
    end function Null_get

!end implementation Null_Field_Expression

!implementation Constant_Field_Expression

    subroutine Constant_set(this, vector)
        class(Constant_Field_Expression), intent(inout) :: this
        real(DP), intent(in) :: vector(:)

        call check_3d_array("Constant_Field_Expression", "vector", vector)
        this%vector = vector
    end subroutine Constant_set

    pure function Constant_get(this, position) result(expression)
        class(Constant_Field_Expression), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: expression(num_dimensions)

        expression = this%vector
    end function Constant_get

!end implementation Constant_Field_Expression

end module class_field_expression
