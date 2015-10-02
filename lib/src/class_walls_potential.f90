module class_walls_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive
use class_periodic_box, only: Abstract_Periodic_Box
use class_floor_penetration, only: Abstract_Floor_Penetration
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Walls_Potential
    private
        real(DP) :: gap
        class(Abstract_Periodic_Box), pointer :: periodic_box
        class(Abstract_Floor_Penetration), pointer :: floor_penetration
    contains
        procedure :: construct => Abstract_Walls_Potential_construct
        procedure :: destroy => Abstract_Walls_Potential_destroy
        procedure :: get_gap => Abstract_Walls_Potential_get_gap
        procedure :: visit => Abstract_Walls_Potential_visit
        procedure, private :: position_from_floor => Abstract_Walls_Potential_position_from_floor
        procedure, private :: position_from_ceiling => &
            Abstract_Walls_Potential_position_from_ceiling
    end type Abstract_Walls_Potential

    type, extends(Abstract_Walls_Potential), public :: Null_Walls_Potential
        contains
        procedure :: construct => Null_Walls_Potential_construct
        procedure :: destroy => Null_Walls_Potential_destroy
        procedure :: get_gap => Null_Walls_Potential_get_gap
        procedure :: visit => Null_Walls_Potential_visit
    end type Null_Walls_Potential

    type, extends(Abstract_Walls_Potential), public :: Concrete_Walls_Potential

    end type Concrete_Walls_Potential

contains

!implementation Abstract_Walls_Potential

    subroutine Abstract_Walls_Potential_construct(this, periodic_box, gap, floor_penetration)
        class(Abstract_Walls_Potential), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: gap
        class(Abstract_Floor_Penetration), target, intent(in) :: floor_penetration


        real(DP) :: box_size(num_dimensions)

        this%periodic_box => periodic_box
        call check_positive("Abstract_Walls_Potential", "gap", gap)
        box_size = this%periodic_box%get_size()
        if (gap > box_size(3)) then
            call error_exit("Abstract_Walls_Potential: gap is too big.")
        end if
        this%gap = gap
        this%floor_penetration => floor_penetration
    end subroutine Abstract_Walls_Potential_construct

    pure function Abstract_Walls_Potential_get_gap(this) result(gap)
        class(Abstract_Walls_Potential), intent(in) :: this
        real(DP) :: gap

        gap = this%gap
    end function Abstract_Walls_Potential_get_gap

    pure subroutine Abstract_Walls_Potential_visit(this, overlap, energy, position, pair_potential)
        class(Abstract_Walls_Potential), intent(in) :: this
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
    end subroutine Abstract_Walls_Potential_visit

    pure function Abstract_Walls_Potential_position_from_floor(this, position) &
        result(position_from_floor)
        class(Abstract_Walls_Potential), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: position_from_floor(num_dimensions)

        position_from_floor = position + [0._DP, 0._DP, this%gap/2._DP]
    end function Abstract_Walls_Potential_position_from_floor

    pure function Abstract_Walls_Potential_position_from_ceiling(this, position) &
        result(position_from_ceiling)
        class(Abstract_Walls_Potential), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: position_from_ceiling(num_dimensions)

        position_from_ceiling = position - [0._DP, 0._DP, this%gap/2._DP]
    end function Abstract_Walls_Potential_position_from_ceiling

    subroutine Abstract_Walls_Potential_destroy(this)
        class(Abstract_Walls_Potential), intent(inout) :: this

        this%floor_penetration => null()
        this%periodic_box => null()
    end subroutine Abstract_Walls_Potential_destroy

!end implementation Abstract_Walls_Potential

!implementation Null_Walls_Potential

    subroutine Null_Walls_Potential_construct(this, periodic_box, gap, floor_penetration)
        class(Null_Walls_Potential), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: gap
        class(Abstract_Floor_Penetration), target, intent(in) :: floor_penetration
    end subroutine Null_Walls_Potential_construct

    pure function Null_Walls_Potential_get_gap(this) result(gap)
        class(Null_Walls_Potential), intent(in) :: this
        real(DP) :: gap
        gap = 0._DP
    end function Null_Walls_Potential_get_gap

    pure subroutine Null_Walls_Potential_visit(this, overlap, energy, position, pair_potential)
        class(Null_Walls_Potential), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        real(DP), intent(in) :: position(:)
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        overlap = .false.
        energy = 0._DP
    end subroutine Null_Walls_Potential_visit

    subroutine Null_Walls_Potential_destroy(this)
        class(Null_Walls_Potential), intent(inout) :: this
    end subroutine Null_Walls_Potential_destroy

!end implementation Null_Walls_Potential

end module class_walls_potential
