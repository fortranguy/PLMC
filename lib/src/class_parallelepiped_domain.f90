module class_parallelepiped_domain

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_periodic_box, only: Abstract_Periodic_Box
use procedures_errors, only: error_exit
use procedures_checks, only: check_3d_array, check_positive
use procedures_parallelepiped_domain, only: point_is_inside

implicit none

private

    type, abstract, public :: Abstract_Parallelepiped_Domain
    private
        real(DP), dimension(num_dimensions) :: origin, size
        class(Abstract_Periodic_Box), pointer :: periodic_box
    contains
        procedure :: construct => Abstract_Parallelepiped_Domain_construct
        procedure, private :: is_boxed => Abstract_Parallelepiped_Domain_is_boxed
        procedure :: destroy => Abstract_Parallelepiped_Domain_destroy
        procedure :: get_volume => Abstract_Parallelepiped_Domain_get_volume
        procedure :: get_vertices => Abstract_Parallelepiped_Domain_get_vertices
        procedure :: is_inside => Abstract_Parallelepiped_Domain_is_inside
        procedure :: overlap => Abstract_Parallelepiped_Domain_overlap
        procedure :: random_position => Abstract_Parallelepiped_Domain_random_position
    end type Abstract_Parallelepiped_Domain

    type, extends(Abstract_Parallelepiped_Domain), public :: Concrete_Parallelepiped_Domain

    end type Concrete_Parallelepiped_Domain

    type, extends(Abstract_Parallelepiped_Domain), public :: Null_Parallelepiped_Domain
    contains
        procedure :: construct => Null_Parallelepiped_Domain_construct
        procedure :: destroy => Null_Parallelepiped_Domain_destroy
        procedure :: get_volume => Null_Parallelepiped_Domain_get_volume
        procedure :: get_vertices => Null_Parallelepiped_Domain_get_vertices
        procedure :: is_inside => Null_Parallelepiped_Domain_is_inside
        procedure :: overlap => Null_Parallelepiped_Domain_overlap
        procedure :: random_position => Null_Parallelepiped_Domain_random_position
    end type Null_Parallelepiped_Domain

    type, extends(Abstract_Parallelepiped_Domain), public :: Box_Parallelepiped_Domain
    contains
        procedure :: construct => Box_Parallelepiped_Domain_construct
    end type Box_Parallelepiped_Domain

contains

!implementation Abstract_Parallelepiped_Domain

    subroutine Abstract_Parallelepiped_Domain_construct(this, periodic_box, origin, size)
        class(Abstract_Parallelepiped_Domain), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: origin(:), size(:)

        this%periodic_box => periodic_box
        call check_3d_array("Abstract_Parallelepiped_Domain", "origin", origin)
        this%origin = this%periodic_box%folded(origin)
        call check_3d_array("Abstract_Parallelepiped_Domain", "size", size)
        call check_positive("Abstract_Parallelepiped_Domain", "size", size)
        this%size = size

        if (.not.this%is_boxed()) then
            call error_exit("Abstract_Parallelepiped_Domain: domain is not boxed.")
        end if
    end subroutine Abstract_Parallelepiped_Domain_construct

    pure logical function Abstract_Parallelepiped_Domain_is_boxed(this) result(is_boxed)
        class(Abstract_Parallelepiped_Domain), intent(in) :: this

        real(DP), dimension(num_dimensions) :: box_origin, corner_minus, corner_plus

        box_origin = 0._DP
        corner_minus = this%origin - this%size/2._DP
        corner_plus = this%origin + this%size/2._DP
        is_boxed = point_is_inside(box_origin, this%periodic_box%get_size(), corner_minus) &
                   .and. point_is_inside(box_origin, this%periodic_box%get_size(), corner_plus)
    end function Abstract_Parallelepiped_Domain_is_boxed


    subroutine Abstract_Parallelepiped_Domain_destroy(this)
        class(Abstract_Parallelepiped_Domain), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Abstract_Parallelepiped_Domain_destroy

    pure real(DP) function Abstract_Parallelepiped_Domain_get_volume(this) result(volume)
        class(Abstract_Parallelepiped_Domain), intent(in) :: this

        volume = product(this%size)
    end function Abstract_Parallelepiped_Domain_get_volume

    pure function Abstract_Parallelepiped_Domain_get_vertices(this, i_vertex) result(vertices)
        real(DP) :: vertices(num_dimensions)
        class(Abstract_Parallelepiped_Domain), intent(in) :: this
        integer, intent(in) :: i_vertex(num_dimensions)

        vertices = this%origin + this%size * (real(i_vertex, DP) - 1.5_DP)
    end function Abstract_Parallelepiped_Domain_get_vertices

    pure logical function Abstract_Parallelepiped_Domain_is_inside(this, position) result(is_inside)
        class(Abstract_Parallelepiped_Domain), intent(in) :: this
        real(DP), intent(in) :: position(:)

        is_inside = point_is_inside(this%origin, this%size, position)
    end function Abstract_Parallelepiped_Domain_is_inside

    pure logical function Abstract_Parallelepiped_Domain_overlap(this, other) result(overlap)
        class(Abstract_Parallelepiped_Domain), intent(in) :: this, other

        logical :: this_vertex_in_other, other_vertex_in_this
        integer :: i_vertex, j_vertex, k_vertex

        overlap = .false.
        do k_vertex = 1, 2
            do j_vertex = 1, 2
                do i_vertex = 1, 2
                    this_vertex_in_other = &
                        other%is_inside(this%get_vertices([i_vertex, j_vertex, k_vertex]))
                    other_vertex_in_this = &
                        this%is_inside(other%get_vertices([i_vertex, j_vertex, k_vertex]))
                    overlap = this_vertex_in_other .or. other_vertex_in_this
                    if (overlap) return
                end do
            end do
        end do
    end function Abstract_Parallelepiped_Domain_overlap

    function Abstract_Parallelepiped_Domain_random_position(this) &
        result(random_position)
        class(Abstract_Parallelepiped_Domain), intent(in) :: this
        real(DP), dimension(num_dimensions) :: random_position

        real(DP) :: rand_3d(num_dimensions)

        call random_number(rand_3d)
        random_position = this%origin + (rand_3d - 0.5_DP) * this%size
    end function Abstract_Parallelepiped_Domain_random_position

!end implementation Abstract_Parallelepiped_Domain

!implementation Null_Parallelepiped_Domain

    subroutine Null_Parallelepiped_Domain_construct(this, periodic_box, origin, size)
        class(Null_Parallelepiped_Domain), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: origin(:), size(:)
    end subroutine Null_Parallelepiped_Domain_construct

    subroutine Null_Parallelepiped_Domain_destroy(this)
        class(Null_Parallelepiped_Domain), intent(inout) :: this
    end subroutine Null_Parallelepiped_Domain_destroy

    pure function Null_Parallelepiped_Domain_get_volume(this) result(volume)
        class(Null_Parallelepiped_Domain), intent(in) :: this
        real(DP) :: volume
        volume = 0._DP
    end function Null_Parallelepiped_Domain_get_volume

    pure function Null_Parallelepiped_Domain_get_vertices(this, i_vertex) result(vertices)
        real(DP) :: vertices(num_dimensions)
        class(Null_Parallelepiped_Domain), intent(in) :: this
        integer, intent(in) :: i_vertex(num_dimensions)
        vertices = 0._DP
    end function Null_Parallelepiped_Domain_get_vertices

    pure function Null_Parallelepiped_Domain_is_inside(this, position) result(is_inside)
        class(Null_Parallelepiped_Domain), intent(in) :: this
        real(DP), intent(in) :: position(:)
        logical :: is_inside
        is_inside = .false.
    end function Null_Parallelepiped_Domain_is_inside

    pure logical function Null_Parallelepiped_Domain_overlap(this, other) result(overlap)
        class(Null_Parallelepiped_Domain), intent(in) :: this
        class(Abstract_Parallelepiped_Domain), intent(in) :: other
        overlap = .false.
    end function Null_Parallelepiped_Domain_overlap

    function Null_Parallelepiped_Domain_random_position(this) &
        result(random_position)
        class(Null_Parallelepiped_Domain), intent(in) :: this
        real(DP), dimension(num_dimensions) :: random_position
        random_position = 0._DP
    end function Null_Parallelepiped_Domain_random_position

!end implementation Null_Parallelepiped_Domain

!implementation Box_Parallelepiped_Domain

    subroutine Box_Parallelepiped_Domain_construct(this, periodic_box, origin, size)
        class(Box_Parallelepiped_Domain), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: origin(:), size(:)

        this%periodic_box => periodic_box
        this%origin = this%periodic_box%folded([0._DP, 0._DP, 0._DP])
        this%size = this%periodic_box%get_size()

        if (.not.this%is_boxed()) then
            call error_exit("Box_Parallelepiped_Domain: domain is not boxed.")
        end if
    end subroutine Box_Parallelepiped_Domain_construct

!end implementation Box_Parallelepiped_Domain

end module class_parallelepiped_domain
