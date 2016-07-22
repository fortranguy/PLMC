module classes_hard_contact

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_dirac_distribution_plus, only: Abstract_Dirac_Distribution_Plus

implicit none

private

    type, abstract, public :: Abstract_Hard_Contact
    private
        class(Abstract_Parallelepiped_Domain), pointer :: accessible_domain => null()
        class(Abstract_Dirac_Distribution_Plus), allocatable :: dirac_plus
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: get_width => Abstract_get_width
        procedure(Abstract_meet), deferred :: meet
        procedure(Abstract_get_beta_pressure_excess), deferred :: get_beta_pressure_excess
    end type Abstract_Hard_Contact

    abstract interface

        pure subroutine Abstract_meet(this, overlap, contact, min_distance, vector)
        import :: DP, Abstract_Hard_Contact
            class(Abstract_Hard_Contact), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: contact
            real(DP), intent(in) :: min_distance !! \( \sigma \)
            real(DP), intent(in) :: vector(:) !! \( \vec{r} \)
        end subroutine Abstract_meet

        pure real(DP) function Abstract_get_beta_pressure_excess(this, contacts)
        import :: DP, Abstract_Hard_Contact
            class(Abstract_Hard_Contact), intent(in) :: this
            real(DP), intent(in) :: contacts
        end function Abstract_get_beta_pressure_excess

    end interface

    type, extends(Abstract_Hard_Contact), public :: XYZ_Hard_Contact
    contains
        procedure :: meet => XYZ_meet
        procedure :: get_beta_pressure_excess => XYZ_get_beta_pressure_excess
    end type XYZ_Hard_Contact

    type, extends(Abstract_Hard_Contact), public :: XY_Hard_Contact
    contains
        procedure :: meet => XY_meet
        procedure :: get_beta_pressure_excess => XY_get_beta_pressure_excess
    end type XY_Hard_Contact

    type, extends(Abstract_Hard_Contact), public :: Null_Hard_Contact
    contains
        procedure :: meet => Null_meet
        procedure :: get_beta_pressure_excess => Null_get_beta_pressure_excess
    end type Null_Hard_Contact

contains

!implementation Abstract_Hard_Contact

    subroutine Abstract_construct(this, accessible_domain, dirac_plus)
        class(Abstract_Hard_Contact), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: accessible_domain
        class(Abstract_Dirac_Distribution_Plus), intent(in) :: dirac_plus

        this%accessible_domain => accessible_domain
        allocate(this%dirac_plus, source=dirac_plus)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Hard_Contact), intent(inout) :: this

        if (allocated(this%dirac_plus)) deallocate(this%dirac_plus)
        this%accessible_domain => null()
    end subroutine Abstract_destroy

    pure real(DP) function Abstract_get_width(this) result(width)
        class(Abstract_Hard_Contact), intent(in) :: this

        width = this%dirac_plus%get_width()
    end function Abstract_get_width

!implementation Abstract_Hard_Contact

!implementation XYZ_Hard_Contact

    !> \[
    !>      \sigma (r = \sigma_+)
    !> \]
    pure subroutine XYZ_meet(this, overlap, contact, min_distance, vector)
        class(XYZ_Hard_Contact), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: contact
        real(DP), intent(in) :: min_distance
        real(DP), intent(in) :: vector(:)

        real(DP) :: distance

        overlap = .false.
        contact = 0._DP
        distance = norm2(vector)
        if (distance < min_distance) then
            overlap = .true.
            return
        end if
        contact = min_distance * this%dirac_plus%get(distance - min_distance)
    end subroutine XYZ_meet

    !> \[
    !>      \frac{1}{3 V} \left\langle \sum_{\mathsf{i} < \mathsf{j}}
    !>          \sigma (r_{\mathsf{i} \mathsf{j}} = \sigma_+) \right\rangle_V
    !> \]
    pure real(DP) function XYZ_get_beta_pressure_excess(this, contacts) result(beta_pressure_excess)
        class(XYZ_Hard_Contact), intent(in) :: this
        real(DP), intent(in) :: contacts

        beta_pressure_excess = 1._DP / (3._DP * product(this%accessible_domain%get_size())) * &
            contacts
    end function XYZ_get_beta_pressure_excess

!end implementation XYZ_Hard_Contact

!implementation XY_Hard_Contact

    !> \[
    !>      \frac{\sigma^2 - z^2}{\sigma} (r = \sigma_+)
    !> \]
    pure subroutine XY_meet(this, overlap, contact, min_distance, vector)
        class(XY_Hard_Contact), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: contact
        real(DP), intent(in) :: min_distance
        real(DP), intent(in) :: vector(:)

        real(DP) :: distance

        overlap = .false.
        contact = 0._DP
        distance = norm2(vector)
        if (distance < min_distance) then
            overlap = .true.
            return
        end if
        contact = (min_distance**2 - vector(3)**2) / min_distance * this%dirac_plus%&
            get(distance - min_distance)
    end subroutine XY_meet

    !> \[
    !>      \frac{1}{2 S H} \left\langle \sum_{\mathsf{i} < \mathsf{j}}
    !>          \frac{\sigma^2 - z_{\mathsf{i} \mathsf{j}}^2}{\sigma}
    !>          (r_{\mathsf{i} \mathsf{j}} = \sigma_+) \right\rangle_{S, H}
    !> \]
    pure real(DP) function XY_get_beta_pressure_excess(this, contacts) result(beta_pressure_excess)
        class(XY_Hard_Contact), intent(in) :: this
        real(DP), intent(in) :: contacts

        beta_pressure_excess = 1._DP / (2._DP * product(this%accessible_domain%get_size())) * &
            contacts
    end function XY_get_beta_pressure_excess

!end implementation XY_Hard_Contact

!implementation Null_Hard_Contact

    pure subroutine Null_meet(this, overlap, contact, min_distance, vector)
        class(Null_Hard_Contact), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: contact
        real(DP), intent(in) :: min_distance
        real(DP), intent(in) :: vector(:)
        overlap = .false.
        contact = 0._DP
    end subroutine Null_meet

    pure real(DP) function Null_get_beta_pressure_excess(this, contacts)&
        result(beta_pressure_excess)
        class(Null_Hard_Contact), intent(in) :: this
        real(DP), intent(in) :: contacts
        beta_pressure_excess = 0._DP
    end function Null_get_beta_pressure_excess


!end implementation Null_Hard_Contact

end module classes_hard_contact
