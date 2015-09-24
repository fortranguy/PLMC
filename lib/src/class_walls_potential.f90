module class_walls_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_floor_penetration, only: Abstract_Floor_Penetration
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Walls_Potential
    private
        real(DP) :: gap
        class(Abstract_Floor_Penetration), pointer :: floor_penetration
        class(Abstract_Pair_Potential), pointer :: pair_potential
    contains
        procedure :: construct => Abstract_Walls_Potential_construct
        procedure :: destroy => Abstract_Walls_Potential_destroy
        procedure :: get_gap => Abstract_Walls_Potential_get_gap
        procedure :: visit => Abstract_Walls_Potential_visit
        procedure, private :: position_from_floor => Abstract_Walls_Potential_position_from_floor
        procedure, private :: position_from_ceiling => &
            Abstract_Walls_Potential_position_from_ceiling
    end type Abstract_Walls_Potential

contains

!implementation Abstract_Walls_Potential

    subroutine Abstract_Walls_Potential_construct(this, gap, floor_penetration, pair_potential)
        class(Abstract_Walls_Potential), intent(out) :: this
        real(DP), intent(in) :: gap
        class(Abstract_Floor_Penetration), target, intent(in) :: floor_penetration
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential

        this%gap = gap
        this%floor_penetration => floor_penetration
        this%pair_potential => pair_potential
    end subroutine Abstract_Walls_Potential_construct

    pure function Abstract_Walls_Potential_get_gap(this) result(gap)
        class(Abstract_Walls_Potential), intent(in) :: this
        real(DP) :: gap

        gap = this%gap
    end function Abstract_Walls_Potential_get_gap

    pure subroutine Abstract_Walls_Potential_visit(this, overlap, energy, position)
        class(Abstract_Walls_Potential), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        real(DP), intent(in) :: position(:)

        real(DP) :: energy_floor, energy_ceiling
        real(DP), dimension(num_dimensions) :: shortest_vector_from_floor, &
            shortest_vector_from_ceiling

        call this%floor_penetration%meet(overlap, shortest_vector_from_floor, &
            this%position_from_floor(position))
        if (overlap) return
        call this%floor_penetration%meet(overlap, shortest_vector_from_ceiling, &
            -this%position_from_ceiling(position))
        if (overlap) return
        call this%pair_potential%meet(norm2(shortest_vector_from_floor), overlap, energy_floor)
        if (overlap) return
        call this%pair_potential%meet(norm2(shortest_vector_from_ceiling), overlap, energy_ceiling)
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

        this%pair_potential => null()
    end subroutine Abstract_Walls_Potential_destroy

!end implementation Abstract_Walls_Potential

end module class_walls_potential
