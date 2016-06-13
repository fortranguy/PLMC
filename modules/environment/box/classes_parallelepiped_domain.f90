module classes_parallelepiped_domain

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_array_size, check_positive
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_visitable_walls, only: Abstract_Visitable_Walls
use procedures_parallelepiped_domain_micro, only: point_is_inside

implicit none

private

    type, abstract, public :: Abstract_Parallelepiped_Domain
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box
    contains
        procedure :: is_boxed => Abstract_is_boxed
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_get_origin), deferred :: get_origin
        procedure(Abstract_get_size), deferred :: get_size
        procedure(Abstract_is_inside), deferred :: is_inside
    end type Abstract_Parallelepiped_Domain

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Parallelepiped_Domain
            class(Abstract_Parallelepiped_Domain), intent(inout) :: this
        end subroutine Abstract_destroy

        pure function Abstract_get_origin(this) result(origin)
        import :: Abstract_Parallelepiped_Domain, DP, num_dimensions
            class(Abstract_Parallelepiped_Domain), intent(in) :: this
            real(DP) :: origin(num_dimensions)
        end function Abstract_get_origin

        pure function Abstract_get_size(this) result(size)
        import :: Abstract_Parallelepiped_Domain, DP, num_dimensions
            class(Abstract_Parallelepiped_Domain), intent(in) :: this
            real(DP) :: size(num_dimensions)
        end function Abstract_get_size

        pure logical function Abstract_is_inside(this, position)
        import :: Abstract_Parallelepiped_Domain, DP
            class(Abstract_Parallelepiped_Domain), intent(in) :: this
            real(DP), intent(in) :: position(:)
        end function Abstract_is_inside

    end interface

    type, extends(Abstract_Parallelepiped_Domain), public :: Concrete_Parallelepiped_Domain
    private
        real(DP), dimension(num_dimensions) :: origin = 0._DP, size = 0._DP
    contains
        procedure :: construct => Parallelepiped_construct
        procedure :: destroy => Parallelepiped_destroy
        procedure :: get_origin => Parallelepiped_get_origin
        procedure :: get_size => Parallelepiped_get_size
        procedure :: is_inside => Parallelepiped_is_inside
    end type Concrete_Parallelepiped_Domain

    type, extends(Abstract_Parallelepiped_Domain), public :: Concrete_Boxed_Domain
    contains
        procedure :: construct => Boxed_construct
        procedure :: destroy => Boxed_destroy
        procedure :: get_origin => Boxed_get_origin
        procedure :: get_size => Boxed_get_size
        procedure :: is_inside => Boxed_is_inside
    end type Concrete_Boxed_Domain

    type, extends(Abstract_Parallelepiped_Domain), public :: Concrete_Walled_Domain
    private
        real(DP) :: size_3
    contains
        procedure :: construct => Walled_construct
        procedure :: destroy => Walled_destroy
        procedure :: get_origin => Walled_get_origin
        procedure :: get_size => Walled_get_size
        procedure :: is_inside => Walled_is_inside
    end type Concrete_Walled_Domain

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

    pure logical function Abstract_is_boxed(this) result(is_boxed)
        class(Abstract_Parallelepiped_Domain), intent(in) :: this

        real(DP), dimension(num_dimensions) :: box_origin, corner_m, corner_p

        box_origin = 0._DP
        corner_m = this%get_origin() - this%get_size()/2._DP
        corner_p = this%get_origin() + this%get_size()/2._DP
        is_boxed = point_is_inside(box_origin, this%periodic_box%get_size(), corner_m) .and. &
            point_is_inside(box_origin, this%periodic_box%get_size(), corner_p)
    end function Abstract_is_boxed

!end implementation Abstract_Parallelepiped_Domain

!implementation Concrete_Parallelepiped_Domain

    subroutine Parallelepiped_construct(this, periodic_box, origin, size)
        class(Concrete_Parallelepiped_Domain), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: origin(:), size(:)

        this%periodic_box => periodic_box
        call check_array_size("Concrete_Parallelepiped_Domain", "origin", origin, num_dimensions)
        this%origin = this%periodic_box%folded(origin)
        call check_array_size("Concrete_Parallelepiped_Domain", "size", size, num_dimensions)
        call check_positive("Concrete_Parallelepiped_Domain", "size", size)
        this%size = size
        if (.not.this%is_boxed()) then
            call error_exit("Parallelepiped_construct: domain is not boxed.")
        end if
    end subroutine Parallelepiped_construct

    subroutine Parallelepiped_destroy(this)
        class(Concrete_Parallelepiped_Domain), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Parallelepiped_destroy

    pure function Parallelepiped_get_origin(this) result(origin)
        class(Concrete_Parallelepiped_Domain), intent(in) :: this
        real(DP) :: origin(num_dimensions)

        origin = this%origin
    end function Parallelepiped_get_origin

    pure function Parallelepiped_get_size(this) result(size)
        class(Concrete_Parallelepiped_Domain), intent(in) :: this
        real(DP) :: size(num_dimensions)

        size = this%size
    end function Parallelepiped_get_size

    pure logical function Parallelepiped_is_inside(this, position) result(is_inside)
        class(Concrete_Parallelepiped_Domain), intent(in) :: this
        real(DP), intent(in) :: position(:)

        is_inside = point_is_inside(this%origin, this%size, position)
    end function Parallelepiped_is_inside

!end implementation Concrete_Parallelepiped_Domain

!implementation Concrete_Boxed_Domain

    subroutine Boxed_construct(this, periodic_box)
        class(Concrete_Boxed_Domain), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box

        this%periodic_box => periodic_box
        if (.not.this%is_boxed()) then
            call error_exit("Boxed_construct: domain is not boxed.")
        end if
    end subroutine Boxed_construct

    subroutine Boxed_destroy(this)
        class(Concrete_Boxed_Domain), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Boxed_destroy

    pure function Boxed_get_origin(this) result(origin)
        class(Concrete_Boxed_Domain), intent(in) :: this
        real(DP) :: origin(num_dimensions)

        origin = 0._DP
    end function Boxed_get_origin

    pure function Boxed_get_size(this) result(size)
        class(Concrete_Boxed_Domain), intent(in) :: this
        real(DP) :: size(num_dimensions)

        size = this%periodic_box%get_size()
    end function Boxed_get_size

    pure logical function Boxed_is_inside(this, position) result(is_inside)
        class(Concrete_Boxed_Domain), intent(in) :: this
        real(DP), intent(in) :: position(:)

        is_inside = .true.
    end function Boxed_is_inside

!end implementation Concrete_Boxed_Domain

!implementation Concrete_Walled_Domain

    subroutine Walled_construct(this, periodic_box, visitable_walls)
        class(Concrete_Walled_Domain), intent(inout) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls

        this%periodic_box => periodic_box
        this%size_3 = visitable_walls%get_max_gap()
        if (.not.this%is_boxed()) then
            call error_exit("Walled_construct: domain is not boxed.")
        end if
    end subroutine Walled_construct

    subroutine Walled_destroy(this)
        class(Concrete_Walled_Domain), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Walled_destroy

    pure function Walled_get_origin(this) result(origin)
        class(Concrete_Walled_Domain), intent(in) :: this
        real(DP) :: origin(num_dimensions)

        origin = 0._DP
    end function Walled_get_origin

    pure function Walled_get_size(this) result(size)
        class(Concrete_Walled_Domain), intent(in) :: this
        real(DP) :: size(num_dimensions)

        size = [reshape(this%periodic_box%get_size(), [2]), this%size_3]
    end function Walled_get_size

    pure logical function Walled_is_inside(this, position) result(is_inside)
        class(Concrete_Walled_Domain), intent(in) :: this
        real(DP), intent(in) :: position(:)

        is_inside = -this%size_3/2._DP <= position(3) .and. position(3) <= this%size_3/2._DP
    end function Walled_is_inside

!end implementation Concrete_Walled_Domain

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
