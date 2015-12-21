module class_potential_expression

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Potential_Expression
    contains
        procedure(Abstract_get), deferred :: get
    end type Abstract_Potential_Expression

    abstract interface

        !> \[u : r \mapsto u(r) \]
        pure real(DP) function Abstract_get(this, distance) result(expression)
        import :: DP, Abstract_Potential_Expression
            class(Abstract_Potential_Expression), intent(in) :: this
            real(DP), intent(in) :: distance
        end function Abstract_get

    end interface

    type, extends(Abstract_Potential_Expression), public :: Null_Potential_Expression
    contains
        procedure :: set => Null_set
        procedure :: get => Null_get
    end type Null_Potential_Expression

    type, extends(Abstract_Potential_Expression), public :: Lennard_Jones_Expression
    private
        real(DP) :: epsilon, sigma
    contains
        procedure :: set => LJ_set
        procedure :: get => LJ_get
    end type Lennard_Jones_Expression

contains

!implementation Lennard_Jones_Expression

    subroutine LJ_set(this, epsilon, sigma)
        class(Lennard_Jones_Expression), intent(inout) :: this
        real(DP), intent(in) :: epsilon, sigma

        call check_positive("Lennard_Jones_Expression: set", "epsilon", epsilon)
        this%epsilon = epsilon
        call check_positive("Lennard_Jones_Expression: set", "sigma", sigma)
        this%sigma =  sigma
    end subroutine LJ_set

    !> \[ u(r) = 4 \epsilon \left[\left(\frac{\sigma}{r}\right)^{12} -
    !>      \left(\frac{\sigma}{r}\right)^6 \right] \]
    pure real(DP) function LJ_get(this, distance) result(expression)
        class(Lennard_Jones_Expression), intent(in) :: this
        real(DP), intent(in) :: distance

        expression = 4._DP * this%epsilon * ((this%sigma/distance)**12 - (this%sigma/distance)**6)
    end function LJ_get

!end implementation Lennard_Jones_Expression

!implementation Null_Potential_Expression

    subroutine Null_set(this)
        class(Null_Potential_Expression), intent(inout) :: this
    end subroutine Null_set

    !> \[ u(r) = 0 \]
    pure real(DP) function Null_get(this, distance) result(expression)
        class(Null_Potential_Expression), intent(in) :: this
        real(DP), intent(in) :: distance
        expression = 0._DP
    end function Null_get

!end implementation Null_Potential_Expression


end module class_potential_expression
