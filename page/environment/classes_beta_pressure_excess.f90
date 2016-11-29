module classes_beta_pressure_excess

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain

implicit none

private

    type, abstract, public :: Abstract_Beta_Pressure_Excess
    private
        class(Abstract_Parallelepiped_Domain), pointer :: accessible_domain => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure(Abstract_get), deferred :: get
    end type Abstract_Beta_Pressure_Excess

    abstract interface

        pure real(DP) function Abstract_get(this, contacts)
        import :: DP, Abstract_Beta_Pressure_Excess
            class(Abstract_Beta_Pressure_Excess), intent(in) :: this
            real(DP), intent(in) :: contacts
        end function Abstract_get

    end interface

    type, extends(Abstract_Beta_Pressure_Excess), public :: XYZ_Beta_Pressure_Excess
    contains
        procedure :: get => XYZ_get
    end type XYZ_Beta_Pressure_Excess

    type, extends(Abstract_Beta_Pressure_Excess), public :: XY_Beta_Pressure_Excess
    contains
        procedure :: get => XY_get
    end type XY_Beta_Pressure_Excess

    type, extends(Abstract_Beta_Pressure_Excess), public :: Null_Beta_Pressure_Excess
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: get => Null_get
    end type Null_Beta_Pressure_Excess

contains

!implementation Abstract_Beta_Pressure_Excess

    subroutine Abstract_construct(this, accessible_domain)
        class(Abstract_Beta_Pressure_Excess), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: accessible_domain

        this%accessible_domain => accessible_domain
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Beta_Pressure_Excess), intent(inout) :: this

        this%accessible_domain => null()
    end subroutine Abstract_destroy

!implementation Abstract_Beta_Pressure_Excess

!implementation XYZ_Beta_Pressure_Excess

    !> \[
    !>      \frac{1}{3 V} \left\langle \sum_{\mathsf{i} < \mathsf{j}}
    !>          \sigma (r_{\mathsf{i} \mathsf{j}} = \sigma_+) \right\rangle_V
    !> \]
    pure real(DP) function XYZ_get(this, contacts) result(beta_pressure_excess)
        class(XYZ_Beta_Pressure_Excess), intent(in) :: this
        real(DP), intent(in) :: contacts

        beta_pressure_excess = 1._DP / (3._DP * product(this%accessible_domain%get_size())) * &
            contacts
    end function XYZ_get

!end implementation XYZ_Beta_Pressure_Excess

!implementation XY_Beta_Pressure_Excess

    !> \[
    !>      \frac{1}{2 S H} \left\langle \sum_{\mathsf{i} < \mathsf{j}}
    !>          \frac{\sigma^2 - z_{\mathsf{i} \mathsf{j}}^2}{\sigma}
    !>          (r_{\mathsf{i} \mathsf{j}} = \sigma_+) \right\rangle_{S, H}
    !> \]
    pure real(DP) function XY_get(this, contacts) result(beta_pressure_excess)
        class(XY_Beta_Pressure_Excess), intent(in) :: this
        real(DP), intent(in) :: contacts

        beta_pressure_excess = 1._DP / (2._DP * product(this%accessible_domain%get_size())) * &
            contacts
    end function XY_get

!end implementation XY_Beta_Pressure_Excess

!implementation Null_Beta_Pressure_Excess

    subroutine Null_construct(this, accessible_domain)
        class(Null_Beta_Pressure_Excess), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: accessible_domain
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Beta_Pressure_Excess), intent(inout) :: this
    end subroutine Null_destroy

    pure real(DP) function Null_get(this, contacts) result(beta_pressure_excess)
        class(Null_Beta_Pressure_Excess), intent(in) :: this
        real(DP), intent(in) :: contacts
        beta_pressure_excess = 0._DP
    end function Null_get

!end implementation Null_Beta_Pressure_Excess

end module classes_beta_pressure_excess
