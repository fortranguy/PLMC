module classes_visitable_walls

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_floor_penetration, only: Abstract_Floor_Penetration
use classes_min_distance, only: Abstract_Min_Distance
use classes_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    !> For the difference between [[Abstract_get_gap_radii]] and [[Abstract_get_gap_centers]],
    !> cf.modules/environment/walls/visitable_walls.tex.
    type, abstract, public :: Abstract_Visitable_Walls
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Floor_Penetration), allocatable :: floor_penetration
        real(DP) :: gap_centers = 0._DP, min_distance = 0._DP
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: are_outside_box => Abstract_are_outside_box
        procedure :: get_gap_radii => Abstract_get_gap_radii
        procedure :: get_gap_centers => Abstract_get_gap_centers !useless?
        procedure :: visit => Abstract_visit
        procedure, private :: position_from_floor => Abstract_position_from_floor
        procedure, private :: position_from_ceiling => Abstract_position_from_ceiling
    end type Abstract_Visitable_Walls

    type, extends(Abstract_Visitable_Walls), public :: Concrete_Visitable_Walls

    end type Concrete_Visitable_Walls

    type, extends(Abstract_Visitable_Walls), public :: Null_Visitable_Walls
        contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: are_outside_box => Null_are_outside_box
        procedure :: get_gap_radii => Null_get_gap
        procedure :: get_gap_centers => Null_get_gap
        procedure :: visit => Null_visit
    end type Null_Visitable_Walls

contains

!implementation Abstract_Visitable_Walls

    subroutine Abstract_construct(this, periodic_box, gap_centers, floor_penetration, min_distance)
        class(Abstract_Visitable_Walls), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: gap_centers
        class(Abstract_Floor_Penetration), intent(in) :: floor_penetration
        class(Abstract_Min_Distance), intent(in) :: min_distance

        this%periodic_box => periodic_box
        call check_positive("Abstract_Visitable_Walls", "gap_centers", gap_centers)
        this%gap_centers = gap_centers
        if (2._DP * floor_penetration%get_height() > this%gap_centers) then
            call error_exit("Abstract_Visitable_Walls: floor_penetration's height is too big.")
        end if
        allocate(this%floor_penetration, source = floor_penetration)
        this%min_distance = min_distance%get()
        if (this%gap_centers < this%min_distance) then
            call error_exit("Abstract_Visitable_Walls overlap: gap_centers is too small.")
        end if
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Visitable_Walls), intent(inout) :: this

        if (allocated(this%floor_penetration)) deallocate(this%floor_penetration)
        this%periodic_box => null()
    end subroutine Abstract_destroy

    pure logical function Abstract_are_outside_box(this) result(are_outside)
        class(Abstract_Visitable_Walls), intent(in) :: this

        real(DP) :: box_size(num_dimensions)

        box_size = this%periodic_box%get_size()
        if (this%get_gap_centers() > box_size(3)) then
            are_outside = .true.
        else
            are_outside = .false.
        end if
    end function Abstract_are_outside_box

    pure real(DP) function Abstract_get_gap_radii(this) result(gap_radii)
        class(Abstract_Visitable_Walls), intent(in) :: this

        gap_radii = this%gap_centers - this%min_distance
    end function Abstract_get_gap_radii

    pure real(DP) function Abstract_get_gap_centers(this) result(gap_centers)
        class(Abstract_Visitable_Walls), intent(in) :: this

        gap_centers = this%gap_centers
    end function Abstract_get_gap_centers

    pure subroutine Abstract_visit(this, overlap, energy, position, pair_potential)
        class(Abstract_Visitable_Walls), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        real(DP), intent(in) :: position(:)
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        real(DP) :: energy_floor, energy_ceiling
        real(DP), dimension(num_dimensions) :: shortest_vector_from_floor, &
            shortest_vector_from_ceiling

        call this%floor_penetration%meet(overlap, shortest_vector_from_floor, &
            this%position_from_floor(position))
        if (overlap) return
        call this%floor_penetration%meet(overlap, shortest_vector_from_ceiling, &
            -this%position_from_ceiling(position))
        if (overlap) return
        call pair_potential%meet(overlap, energy_floor, norm2(shortest_vector_from_floor))
        if (overlap) return
        call pair_potential%meet(overlap, energy_ceiling, norm2(shortest_vector_from_ceiling))
        if (overlap) return
        energy = energy_floor + energy_ceiling
    end subroutine Abstract_visit

    pure function Abstract_position_from_floor(this, position) result(position_from_floor)
        class(Abstract_Visitable_Walls), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: position_from_floor(num_dimensions)

        position_from_floor = position + [0._DP, 0._DP, this%gap_centers/2._DP]
    end function Abstract_position_from_floor

    pure function Abstract_position_from_ceiling(this, position) result(position_from_ceiling)
        class(Abstract_Visitable_Walls), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: position_from_ceiling(num_dimensions)

        position_from_ceiling = position - [0._DP, 0._DP, this%gap_centers/2._DP]
    end function Abstract_position_from_ceiling

!end implementation Abstract_Visitable_Walls

!implementation Null_Visitable_Walls

    subroutine Null_construct(this, periodic_box, gap_centers, floor_penetration, min_distance)
        class(Null_Visitable_Walls), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: gap_centers
        class(Abstract_Floor_Penetration), intent(in) :: floor_penetration
        class(Abstract_Min_Distance), intent(in) :: min_distance
    end subroutine Null_construct

    pure logical function Null_are_outside_box(this) result(are_outside)
        class(Null_Visitable_Walls), intent(in) :: this
        are_outside = .false.
    end function Null_are_outside_box

    pure real(DP) function Null_get_gap(this) result(gap)
        class(Null_Visitable_Walls), intent(in) :: this
        gap = 0._DP
    end function Null_get_gap

    pure subroutine Null_visit(this, overlap, energy, position, pair_potential)
        class(Null_Visitable_Walls), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        real(DP), intent(in) :: position(:)
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        overlap = .false.
        energy = 0._DP
    end subroutine Null_visit

    subroutine Null_destroy(this)
        class(Null_Visitable_Walls), intent(inout) :: this
    end subroutine Null_destroy

!end implementation Null_Visitable_Walls

end module classes_visitable_walls
