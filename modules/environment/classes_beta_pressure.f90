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

contains

    subroutine Abstract_set(this, density, excess)
        class(Abstract_Beta_Pressure), intent(inout) :: this
        real(DP), intent(in) :: density !! if with walls: \( \frac{N}{H S} \)
        real(DP), intent(in) :: excess

        call check_positive("Abstract_Beta_Pressure: set", "density", density)
        this%density = density !smarter?
        this%excess = excess !what check?
    end subroutine Abstract_set

    !> \( \beta p_\text = \rho + \beta p_\text{ex} \)
    pure real(DP) function Abstract_get(this) result(beta_pressure)
        class(Abstract_Beta_Pressure), intent(in) :: this

        beta_pressure = this%density + this%excess
    end function Abstract_get

end module classes_beta_pressure
