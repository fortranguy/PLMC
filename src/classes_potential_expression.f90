module classes_potential_expression

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Potential_Expression
    contains
        procedure(Abstract_get), deferred :: get
    end type Abstract_Potential_Expression

    abstract interface

        !> \[r \mapsto u(r) \]
        pure real(DP) function Abstract_get(this, distance) result(expression)
        import :: DP, Abstract_Potential_Expression
            class(Abstract_Potential_Expression), intent(in) :: this
            real(DP), intent(in) :: distance !! \( r \)
        end function Abstract_get

    end interface

    type, extends(Abstract_Potential_Expression), public :: Null_Potential_Expression
    contains
        procedure :: get => Null_get
    end type Null_Potential_Expression

    type, extends(Abstract_Potential_Expression), public :: Lennard_Jones_Expression
    private
        real(DP) :: epsilon = 0._DP
        real(DP) :: sigma = 0._DP
    contains
        procedure :: set => LJ_set
        procedure :: get => LJ_get
    end type Lennard_Jones_Expression

    type, extends(Abstract_Potential_Expression), public :: Yukawa_Expression
    private
        real(DP) :: epsilon = 0._DP
        real(DP) :: alpha = 0._DP
    contains
        procedure :: set => Yukawa_set
        procedure :: get => Yukawa_get
    end type Yukawa_Expression

contains

!implementation Lennard_Jones_Expression

    subroutine LJ_set(this, epsilon, sigma)
        class(Lennard_Jones_Expression), intent(inout) :: this
        real(DP), intent(in) :: epsilon !! \( \epsilon \)
        real(DP), intent(in) :: sigma !! \( \sigma \)

        call check_positive("Lennard_Jones_Expression: set", "epsilon", epsilon)
        this%epsilon = epsilon
        call check_positive("Lennard_Jones_Expression: set", "sigma", sigma)
        this%sigma =  sigma
    end subroutine LJ_set

    !> \[
    !>      r \mapsto 4 \epsilon \left[\left(\frac{\sigma}{r}\right)^{12} -
    !>          \left(\frac{\sigma}{r}\right)^6 \right]
    !> \]
    pure real(DP) function LJ_get(this, distance) result(expression)
        class(Lennard_Jones_Expression), intent(in) :: this
        real(DP), intent(in) :: distance

        expression = 4._DP * this%epsilon * ((this%sigma/distance)**12 - (this%sigma/distance)**6)
    end function LJ_get

!end implementation Lennard_Jones_Expression

!implementation Yukawa_Expression

    subroutine Yukawa_set(this, epsilon, alpha)
        class(Yukawa_Expression), intent(inout) :: this
        real(DP), intent(in) :: epsilon !! \( \epsilon \)
        real(DP), intent(in) :: alpha !! \( \alpha \)

        this%epsilon = epsilon
        call check_positive("Yukawa_Expression: set", "alpha", alpha)
        this%alpha = alpha
    end subroutine Yukawa_set

    !> \[ r \mapsto \frac{\epsilon}{r} \exp(-\alpha r) \]
    pure real(DP) function Yukawa_get(this, distance) result(expression)
        class(Yukawa_Expression), intent(in) :: this
        real(DP), intent(in) :: distance

        expression = this%epsilon / distance * exp(-this%alpha * distance)
    end function Yukawa_get

!end implementation Yukawa_Expression

!implementation Null_Potential_Expression

    !> \[ r \mapto 0 \]
    pure real(DP) function Null_get(this, distance) result(expression)
        class(Null_Potential_Expression), intent(in) :: this
        real(DP), intent(in) :: distance
        expression = 0._DP
    end function Null_get

!end implementation Null_Potential_Expression


end module classes_potential_expression
