module class_potential_expression

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use types_potential_parameters, only: Abstract_Potential_Parameters, Lennard_Jones_Parameters

implicit none

private

    type, abstract, public :: Abstract_Potential_Expression
    contains
        procedure(Abstract_Potential_Expression_set), deferred :: set
        procedure(Abstract_Potential_Expression_get), deferred :: get
    end type Abstract_Potential_Expression

    abstract interface

        subroutine Abstract_Potential_Expression_set(this, parameters)
        import :: Abstract_Potential_Parameters, Abstract_Potential_Expression
            class(Abstract_Potential_Expression), intent(inout) :: this
            class(Abstract_Potential_Parameters), intent(in) :: parameters
        end subroutine Abstract_Potential_Expression_set

        pure function Abstract_Potential_Expression_get(this, distance) result(field_expression)
        import :: DP, Abstract_Potential_Expression
            class(Abstract_Potential_Expression), intent(in) :: this
            real(DP), intent(in) :: distance
            real(DP) :: field_expression
        end function Abstract_Potential_Expression_get

    end interface

    type, extends(Abstract_Potential_Expression), public :: Null_Potential_Expression
    contains
        procedure :: set => Null_Potential_Expression_set
        procedure :: get => Null_Potential_Expression_get
    end type Null_Potential_Expression

    type, extends(Abstract_Potential_Expression), public :: Lennard_Jones_Expression
    private
        type(Lennard_Jones_Parameters) :: parameters
    contains
        procedure :: set => Lennard_Jones_Expression_set
        procedure :: get => Lennard_Jones_Expression_get
    end type Lennard_Jones_Expression

contains

!implementation Null_Potential_Expression

    subroutine Null_Potential_Expression_set(this, parameters)
        class(Null_Potential_Expression), intent(inout) :: this
        class(Abstract_Potential_Parameters), intent(in) :: parameters
    end subroutine Null_Potential_Expression_set

    pure function Null_Potential_Expression_get(this, distance) result(expression)
        class(Null_Potential_Expression), intent(in) :: this
        real(DP), intent(in) :: distance
        real(DP) :: expression
        expression = 0._DP
    end function Null_Potential_Expression_get

!end implementation Null_Potential_Expression

!implemetation Lennard_Jones_Expression

    subroutine Lennard_Jones_Expression_set(this, parameters)
        class(Lennard_Jones_Expression), intent(inout) :: this
        class(Abstract_Potential_Parameters), intent(in) :: parameters

        select type (parameters)
            type is (Lennard_Jones_Parameters)
                this%parameters%epsilon = parameters%epsilon
                this%parameters%sigma =  parameters%sigma
            class default
                call error_exit("Lennard_Jones_Expression: no parameters were given.")
        end select
    end subroutine Lennard_Jones_Expression_set

    pure function Lennard_Jones_Expression_get(this, distance) result(expression)
        class(Lennard_Jones_Expression), intent(in) :: this
        real(DP), intent(in) :: distance
        real(DP) :: expression

        real(DP) :: epsilon, sigma

        epsilon = this%parameters%epsilon
        sigma = this%parameters%sigma
        expression = 4._DP * epsilon * ((sigma/distance)**12 - (sigma/distance)**6)
    end function Lennard_Jones_Expression_get

!end implemetation Lennard_Jones_Expression

end module class_potential_expression
