module class_potential_expression

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Potential_Expression
    contains
        procedure(Abstract_Potential_Expression_get), deferred :: get
    end type Abstract_Potential_Expression

    abstract interface

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
        real(DP) :: epsilon, sigma
    contains
        procedure :: set => Lennard_Jones_Expression_set
        procedure :: get => Lennard_Jones_Expression_get
    end type Lennard_Jones_Expression

contains

!implementation Null_Potential_Expression

    subroutine Null_Potential_Expression_set(this)
        class(Null_Potential_Expression), intent(inout) :: this
    end subroutine Null_Potential_Expression_set

    pure function Null_Potential_Expression_get(this, distance) result(expression)
        class(Null_Potential_Expression), intent(in) :: this
        real(DP), intent(in) :: distance
        real(DP) :: expression
        expression = 0._DP
    end function Null_Potential_Expression_get

!end implementation Null_Potential_Expression

!implemetation Lennard_Jones_Expression

    subroutine Lennard_Jones_Expression_set(this, epsilon, sigma)
        class(Lennard_Jones_Expression), intent(inout) :: this
        real(DP), intent(in) :: epsilon, sigma

        call check_positive("Lennard_Jones_Expression", "epsilon", epsilon)
        this%epsilon = epsilon
        call check_positive("Lennard_Jones_Expression", "sigma", sigma)
        this%sigma =  sigma
    end subroutine Lennard_Jones_Expression_set

    pure function Lennard_Jones_Expression_get(this, distance) result(expression)
        class(Lennard_Jones_Expression), intent(in) :: this
        real(DP), intent(in) :: distance
        real(DP) :: expression

        expression = 4._DP * this%epsilon * ((this%sigma/distance)**12 - (this%sigma/distance)**6)
    end function Lennard_Jones_Expression_get

!end implemetation Lennard_Jones_Expression

end module class_potential_expression
