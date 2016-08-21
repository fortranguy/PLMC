module classes_beta_pressure

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Beta_Pressure
    private
        real(DP) :: density
        real(DP) :: excess
    contains
        procedure :: set => Abstract_set
        procedure :: get => Abstract_get
    end type Abstract_Beta_Pressure

    type, extends(Abstract_Beta_Pressure), public :: Concrete_Beta_Pressure

    end type Concrete_Beta_Pressure

    type, extends(Abstract_Beta_Pressure), public :: Null_Beta_Pressure
    contains
        procedure :: set => Null_set
        procedure :: get => Null_get
    end type Null_Beta_Pressure

contains

!implementation Abstract_Beta_Pressure

    subroutine Abstract_set(this, density, excess)
        class(Abstract_Beta_Pressure), intent(inout) :: this
        real(DP), intent(in) :: density !! if with walls: \( \frac{N}{H S} \)
        real(DP), intent(in) :: excess

        call check_positive("Abstract_Beta_Pressure: set", "density", density)
        this%density = density !smarter?
        this%excess = excess !what check?
    end subroutine Abstract_set

    !> \( \beta p = \rho + \beta p_\text{ex} \)
    pure real(DP) function Abstract_get(this) result(beta_pressure)
        class(Abstract_Beta_Pressure), intent(in) :: this

        beta_pressure = this%density + this%excess
    end function Abstract_get

!implementation Abstract_Beta_Pressure

!implementation Null_Beta_Pressure

    subroutine Null_set(this, density, excess)
        class(Null_Beta_Pressure), intent(inout) :: this
        real(DP), intent(in) :: density
        real(DP), intent(in) :: excess
    end subroutine Null_set

    pure real(DP) function Null_get(this) result(beta_pressure)
        class(Null_Beta_Pressure), intent(in) :: this
        beta_pressure = 0._DP
    end function Null_get

!implementation Null_Beta_Pressure

end module classes_beta_pressure
