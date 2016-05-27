module classes_random_orientation

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_normal_random_number, only: normal_random_number
use classes_random_coordinates, only: Abstract_Random_Coordinates

implicit none

private

    type, extends(Abstract_Random_Coordinates), public :: Concrete_Random_Orientation
    private
        logical, allocatable :: have_orientations(:)
    contains
        procedure :: construct => Concrete_construct
        procedure :: destroy => Concrete_destroy
        procedure :: get => Concrete_get
    end type Concrete_Random_Orientation

contains

    subroutine Concrete_construct(this, have_orientations)
        class(Concrete_Random_Orientation), intent(out) :: this
        logical, intent(in) :: have_orientations(:)

        allocate(this%have_orientations, source=have_orientations)
    end subroutine Concrete_construct

    subroutine Concrete_destroy(this)
        class(Concrete_Random_Orientation), intent(inout) :: this

        if (allocated(this%have_orientations)) deallocate(this%have_orientations)
    end subroutine Concrete_destroy

    !> From SMAC, Algorithm 1.23, p. 43
    function Concrete_get(this, i_component) result(random_orientation)
        class(Concrete_Random_Orientation), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP) :: random_orientation(num_dimensions)

        if (this%have_orientations(i_component)) then
            call normal_random_number(random_orientation)
            random_orientation = random_orientation / norm2(random_orientation)
        else
            random_orientation = 0._DP
        end if
    end function Concrete_get

end module classes_random_orientation
