module classes_parallelepiped_domain

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_errors, only: warning_continue
use procedures_checks, only: check_array_size, check_positive

implicit none

private

    type, abstract, public :: Abstract_Parallelepiped_Domain
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box
        real(DP), dimension(num_dimensions) :: origin = 0._DP, size = 0._DP
    contains
        procedure :: construct => Abstract_construct
        procedure, private :: is_boxed => Abstract_is_boxed
        procedure :: destroy => Abstract_destroy
        procedure :: get_origin => Abstract_get_origin
        procedure :: get_size => Abstract_get_size
        procedure :: is_inside => Abstract_is_inside
    end type Abstract_Parallelepiped_Domain

    type, extends(Abstract_Parallelepiped_Domain), public :: Concrete_Parallelepiped_Domain

    end type Concrete_Parallelepiped_Domain

    type, extends(Abstract_Parallelepiped_Domain), public :: Concrete_Box_Domain
    contains
        procedure :: construct => Box_construct
        procedure :: get_origin => Box_get_origin
        procedure :: get_size => Box_get_size
        procedure :: is_inside => Box_is_inside
    end type Concrete_Box_Domain

    type, extends(Abstract_Parallelepiped_Domain), public :: Null_Parallelepiped_Domain
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: get_origin => Null_get_origin
        procedure :: get_size => Null_get_size
        procedure :: is_inside => Null_is_inside
    end type Null_Parallelepiped_Domain

contains

!implementation Abstract_Parallelepiped_Domain

    subroutine Abstract_construct(this, periodic_box, origin, size)
        class(Abstract_Parallelepiped_Domain), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: origin(:), size(:)

        this%periodic_box => periodic_box
        call check_array_size("Abstract_Parallelepiped_Domain", "origin", origin, num_dimensions)
        this%origin = this%periodic_box%folded(origin)
        call check_array_size("Abstract_Parallelepiped_Domain", "size", size, num_dimensions)
        call check_positive("Abstract_Parallelepiped_Domain", "size", size)
        this%size = size

        if (.not.this%is_boxed()) then
            call warning_continue("Abstract_Parallelepiped_Domain: domain is not boxed.")
        end if
    end subroutine Abstract_construct

    pure logical function Abstract_is_boxed(this) result(is_boxed)
        class(Abstract_Parallelepiped_Domain), intent(in) :: this

        real(DP), dimension(num_dimensions) :: box_origin, corner_minus, corner_plus

        box_origin = 0._DP
        corner_minus = this%origin - this%size/2._DP
        corner_plus = this%origin + this%size/2._DP
        is_boxed = point_is_inside(box_origin, this%periodic_box%get_size(), corner_minus) &
                   .and. point_is_inside(box_origin, this%periodic_box%get_size(), corner_plus)
    end function Abstract_is_boxed

    subroutine Abstract_destroy(this)
        class(Abstract_Parallelepiped_Domain), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Abstract_destroy

    pure function Abstract_get_origin(this) result(origin)
        class(Abstract_Parallelepiped_Domain), intent(in) :: this
        real(DP) :: origin(num_dimensions)

        origin = this%origin
    end function Abstract_get_origin

    pure function Abstract_get_size(this) result(size)
        class(Abstract_Parallelepiped_Domain), intent(in) :: this
        real(DP) :: size(num_dimensions)

        size = this%size
    end function Abstract_get_size

    pure logical function Abstract_is_inside(this, position) result(is_inside)
        class(Abstract_Parallelepiped_Domain), intent(in) :: this
        real(DP), intent(in) :: position(:)

        is_inside = point_is_inside(this%origin, this%size, position)
    end function Abstract_is_inside

    pure logical function point_is_inside(box_origin, box_size, point)
        real(DP), intent(in) :: box_origin(:), box_size(:), point(:)

        real(DP), dimension(num_dimensions) :: box_corner, point_from_corner

        box_corner = box_origin - box_size/2._DP
        point_from_corner = point - box_corner
        point_is_inside = all(0._DP <= point_from_corner .and. point_from_corner <= box_size)
    end function point_is_inside

!end implementation Abstract_Parallelepiped_Domain

!implementation Concrete_Box_Domain

    subroutine Box_construct(this, periodic_box, origin, size)
        class(Concrete_Box_Domain), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: origin(:), size(:)

        this%periodic_box => periodic_box
    end subroutine Box_construct

    pure function Box_get_origin(this) result(origin)
        class(Concrete_Box_Domain), intent(in) :: this
        real(DP) :: origin(num_dimensions)

        origin = 0._DP
    end function Box_get_origin

    pure function Box_get_size(this) result(size)
        class(Concrete_Box_Domain), intent(in) :: this
        real(DP) :: size(num_dimensions)

        size = this%periodic_box%get_size()
    end function Box_get_size

    pure logical function Box_is_inside(this, position) result(is_inside)
        class(Concrete_Box_Domain), intent(in) :: this
        real(DP), intent(in) :: position(:)

        is_inside = .true.
    end function Box_is_inside

!end implementation Concrete_Box_Domain

!implementation Null_Parallelepiped_Domain

    subroutine Null_construct(this, periodic_box, origin, size)
        class(Null_Parallelepiped_Domain), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: origin(:), size(:)
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Parallelepiped_Domain), intent(inout) :: this
    end subroutine Null_destroy

    pure function Null_get_origin(this) result(origin)
        class(Null_Parallelepiped_Domain), intent(in) :: this
        real(DP) :: origin(num_dimensions)
        origin = 0._DP
    end function Null_get_origin

    pure function Null_get_size(this) result(size)
        class(Null_Parallelepiped_Domain), intent(in) :: this
        real(DP) :: size(num_dimensions)
        size = 0._DP
    end function Null_get_size

    pure function Null_is_inside(this, position) result(is_inside)
        class(Null_Parallelepiped_Domain), intent(in) :: this
        real(DP), intent(in) :: position(:)
        logical :: is_inside
        is_inside = .false.
    end function Null_is_inside

!end implementation Null_Parallelepiped_Domain

end module classes_parallelepiped_domain
