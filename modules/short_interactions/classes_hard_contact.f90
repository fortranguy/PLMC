module classes_hard_contact

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_half_distribution, only: Abstract_Half_Distribution

implicit none

private

    type, abstract, public :: Abstract_Hard_Contact
    private
        class(Abstract_Half_Distribution), allocatable :: half_distribution
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure(Abstract_meet), deferred :: meet
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

    end interface

    type, extends(Abstract_Hard_Contact), public :: XYZ_Hard_Contact
    contains
        procedure :: meet => XYZ_meet
    end type XYZ_Hard_Contact

    type, extends(Abstract_Hard_Contact), public :: XY_Hard_Contact
    contains
        procedure :: meet => XY_meet
    end type XY_Hard_Contact

contains

!implementation Abstract_Hard_Contact

    subroutine Abstract_construct(this, half_distribution)
        class(Abstract_Hard_Contact), intent(out) :: this
        class(Abstract_Half_Distribution), intent(in) :: half_distribution

        allocate(this%half_distribution, source=half_distribution)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Hard_Contact), intent(inout) :: this

        if (allocated(this%half_distribution)) deallocate(this%half_distribution)
    end subroutine Abstract_destroy

!implementation Abstract_Hard_Contact

!implementation XYZ_Hard_Contact

    !> \[
    !>      \sigma \lvert r = \sigma)
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
        contact = min_distance * this%half_distribution%get(distance - min_distance)
    end subroutine XYZ_meet

!end implementation XYZ_Hard_Contact

!implementation XY_Hard_Contact

    !> \[
    !>      \frac{\sigma^2 - z^2}{\sigma} \lvert r = \sigma)
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
        contact = (min_distance**2 - vector(3)**2) / min_distance * this%half_distribution%&
            get(distance - min_distance)
    end subroutine XY_meet

!end implementation XY_Hard_Contact

end module classes_hard_contact
