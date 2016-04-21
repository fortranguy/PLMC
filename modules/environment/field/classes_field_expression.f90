module classes_field_expression

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use procedures_errors, only: error_exit
use procedures_checks, only: check_array_size, check_positive
use classes_permittivity, only: Abstract_Permittivity
use procedures_field_expression_micro, only: plate_expression

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

    type, extends(Abstract_Field_Expression), public :: Constant_Field_Expression
    private
        real(DP) :: vector(num_dimensions) = 0._DP
    contains
        procedure :: set => Constant_set
        procedure :: get => Constant_get
    end type Constant_Field_Expression

    !> This field is produced by 2 charged plates of opposite charge.
    type, extends(Abstract_Field_Expression), public :: Centered_Plates_Expression
    private
        real(DP) :: coulomb = 0._DP
        real(DP) :: size_x = 0._DP
        real(DP), dimension(2) :: center_lower = 0._DP, center_upper = 0._DP
        real(DP) :: surface_density = 0._DP
    contains
        procedure :: set => Plates_set
        procedure :: get => Plates_get
    end type Centered_Plates_Expression

    type, extends(Abstract_Field_Expression), public :: Null_Field_Expression
    contains
        procedure :: set => Null_set
        procedure :: get => Null_get
    end type Null_Field_Expression

contains

!implementation Constant_Field_Expression

    subroutine Constant_set(this, vector)
        class(Constant_Field_Expression), intent(inout) :: this
        real(DP), intent(in) :: vector(:)

        call check_array_size("Constant_Field_Expression", "vector", vector, num_dimensions)
        this%vector = vector
    end subroutine Constant_set

    pure function Constant_get(this, position) result(expression)
        class(Constant_Field_Expression), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: expression(num_dimensions)

        expression = this%vector
    end function Constant_get

!end implementation Constant_Field_Expression

!implementation Centered_Plates_Expression

    subroutine Plates_set(this, permittivity, gap, size_x, surface_density)
        class(Centered_Plates_Expression), intent(inout) :: this
        class(Abstract_Permittivity), intent(in) :: permittivity
        real(DP), intent(in) :: gap, size_x
        real(DP), intent(in) :: surface_density

        this%coulomb = 1._DP / (4._DP*PI * permittivity%get())
        this%surface_density = surface_density
        call check_positive("Centered_Plates_Expression: set", "size_x", size_x)
        this%size_x = size_x
        call check_positive("Centered_Plates_Expression: set", "gap", gap)
        this%center_lower = [0._DP, -gap/2._DP]
        this%center_upper = [0._DP, +gap/2._DP]
    end subroutine Plates_set

    !> Let an infinitely thin plate of charge density \( \sigma \)
    !> which is finite in \( \vec{e}_x \) direction but not in \( \vec{e}_y \),
    !> cf. modules/environment/plate_field.tex.
    !> This plate creates an electric field at \( (x, z) \) of expression:
    !> \[
    !>      \vec{E}(x, z) = \frac{1}{4\pi \epsilon_0} \sigma \vec{b}(x, z)
    !> \]
    !> where \( \vec{b}(x, z) \) is defined in
    !> [[procedures_field_expression_micro:plate_expression]].
    pure function Plates_get(this, position) result(expression)
        class(Centered_Plates_Expression), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: expression(num_dimensions)

        real(DP), dimension(2) :: position_lower_13, position_upper_13

        position_lower_13 = [position(1), position(3)] - this%center_lower
        position_upper_13 = [position(1), position(3)] - this%center_upper
        expression = this%coulomb * this%surface_density * &
            (plate_expression(this%size_x, position_lower_13) - &
             plate_expression(this%size_x, position_upper_13))
    end function Plates_get

!end implementation Centered_Plates_Expression

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

end module classes_field_expression
