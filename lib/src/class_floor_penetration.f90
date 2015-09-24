module class_floor_penetration

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions

implicit none

private

    type, abstract, public :: Abstract_Floor_Penetration
    contains
        procedure :: meet => Abstract_Floor_Penetration_meet
    end type Abstract_Floor_Penetration

contains

    pure subroutine Abstract_Floor_Penetration_meet(this, overlap, shortest_vector_from_floor, &
        position_from_floor)
        class(Abstract_Floor_Penetration), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: shortest_vector_from_floor(num_dimensions)
        real(DP), intent(in) :: position_from_floor(num_dimensions)

        shortest_vector_from_floor = [0._DP, 0._DP, position_from_floor(3)]
        if (shortest_vector_from_floor(3) < 0._DP) then
            overlap = .true.
        else
            overlap = .false.
        end if
    end subroutine Abstract_Floor_Penetration_meet

end module class_floor_penetration
