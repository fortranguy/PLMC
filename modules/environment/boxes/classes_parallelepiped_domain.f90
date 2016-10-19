module classes_parallelepiped_domain

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_array_size, check_positive
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_visitable_walls, only: Abstract_Visitable_Walls
use procedures_parallelepiped_domain_micro, only: point_is_inside_box, is_inside_included, &
    is_inside_excluded

implicit none

private

    type, abstract, public :: Abstract_Parallelepiped_Domain
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_get_origin), deferred :: get_origin
        procedure(Abstract_get_size), deferred :: get_size
        procedure :: is_boxed => Abstract_is_boxed
        procedure :: is_inside => Abstract_is_inside
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

    end interface

    type, extends(Abstract_Parallelepiped_Domain), public :: Concrete_Parallelepiped_Domain
    private
        real(DP), dimension(num_dimensions) :: origin = 0._DP, size = 0._DP
    contains
        procedure :: construct => Parallelepiped_construct
        procedure :: destroy => Parallelepiped_destroy
        procedure :: get_origin => Parallelepiped_get_origin
        procedure :: get_size => Parallelepiped_get_size
    end type Concrete_Parallelepiped_Domain

    type, extends(Abstract_Parallelepiped_Domain), public :: Boxed_Parallelepiped_Domain
    contains
        procedure :: construct => Boxed_construct
        procedure :: destroy => Boxed_destroy
        procedure :: get_origin => Boxed_get_origin
        procedure :: get_size => Boxed_get_size
    end type Boxed_Parallelepiped_Domain

    type, extends(Abstract_Parallelepiped_Domain), public :: Walled_Parallelepiped_Domain
    private
        real(DP) :: size_3
    contains
        procedure :: construct => Walled_construct
        procedure :: destroy => Walled_destroy
        procedure :: get_origin => Walled_get_origin
        procedure :: get_size => Walled_get_size
    end type Walled_Parallelepiped_Domain

    type, extends(Abstract_Parallelepiped_Domain), public :: Null_Parallelepiped_Domain
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: get_origin => Null_get_origin
        procedure :: get_size => Null_get_size
        procedure :: is_boxed => Null_is_boxed
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
        is_boxed = point_is_inside_box(box_origin, this%periodic_box%get_size(), corner_m, &
            is_inside_included) .and. point_is_inside_box(box_origin, this%periodic_box%get_size(),&
            corner_p, is_inside_included)
    end function Abstract_is_boxed

    pure logical function Abstract_is_inside(this, position) result(is_inside)
        class(Abstract_Parallelepiped_Domain), intent(in) :: this
        real(DP), intent(in) :: position(:)

        is_inside = point_is_inside_box(this%get_origin(), this%get_size(), position, &
            is_inside_excluded)
    end function Abstract_is_inside

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

!end implementation Concrete_Parallelepiped_Domain

!implementation Boxed_Parallelepiped_Domain

    subroutine Boxed_construct(this, periodic_box)
        class(Boxed_Parallelepiped_Domain), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box

        this%periodic_box => periodic_box
    end subroutine Boxed_construct

    subroutine Boxed_destroy(this)
        class(Boxed_Parallelepiped_Domain), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Boxed_destroy

    pure function Boxed_get_origin(this) result(origin)
        class(Boxed_Parallelepiped_Domain), intent(in) :: this
        real(DP) :: origin(num_dimensions)

        origin = 0._DP
    end function Boxed_get_origin

    pure function Boxed_get_size(this) result(size)
        class(Boxed_Parallelepiped_Domain), intent(in) :: this
        real(DP) :: size(num_dimensions)

        size = this%periodic_box%get_size()
    end function Boxed_get_size

!end implementation Boxed_Parallelepiped_Domain

!implementation Walled_Parallelepiped_Domain

    subroutine Walled_construct(this, periodic_box, visitable_walls)
        class(Walled_Parallelepiped_Domain), intent(inout) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls

        this%periodic_box => periodic_box
        this%size_3 = visitable_walls%get_gap_radii()
    end subroutine Walled_construct

    subroutine Walled_destroy(this)
        class(Walled_Parallelepiped_Domain), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Walled_destroy

    pure function Walled_get_origin(this) result(origin)
        class(Walled_Parallelepiped_Domain), intent(in) :: this
        real(DP) :: origin(num_dimensions)

        origin = 0._DP
    end function Walled_get_origin

    pure function Walled_get_size(this) result(size)
        class(Walled_Parallelepiped_Domain), intent(in) :: this
        real(DP) :: size(num_dimensions)

        size = [reshape(this%periodic_box%get_size(), [2]), this%size_3]
    end function Walled_get_size

!end implementation Walled_Parallelepiped_Domain

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

    pure logical function Null_is_boxed(this) result(is_boxed)
        class(Null_Parallelepiped_Domain), intent(in) :: this
        is_boxed = .true.
    end function Null_is_boxed

    pure function Null_is_inside(this, position) result(is_inside)
        class(Null_Parallelepiped_Domain), intent(in) :: this
        real(DP), intent(in) :: position(:)
        logical :: is_inside
        is_inside = .false.
    end function Null_is_inside

!end implementation Null_Parallelepiped_Domain

end module classes_parallelepiped_domain
