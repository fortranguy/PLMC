module procedures_coordinates_copier_factory

use procedures_errors, only: error_exit
use classes_random_coordinates, only: Abstract_Random_Coordinates
use classes_coordinates_copier, only: Abstract_Coordinates_Copier, &
    Random_Filling_Coordinates_Copier, Null_Coordinates_Copier

implicit none

private
public :: create_position, create_orientation, destroy

contains

    subroutine create_position(position_copier, random_position, have_positions)
        class(Abstract_Coordinates_Copier), allocatable, intent(out) :: position_copier
        class(Abstract_Random_Coordinates), intent(in) :: random_position
        logical, intent(in) :: have_positions(:)

        if (any(have_positions)) then
            allocate(Random_Filling_Coordinates_Copier :: position_copier)
        else
            allocate(Null_Coordinates_Copier :: position_copier)
        end if
        select type (position_copier)
            type is (Random_Filling_Coordinates_Copier)
                call position_copier%construct(random_position, have_positions)
            type is (Null_Coordinates_Copier)
            class default
                call error_exit("procedures_coordinates_copier_factory: create_position: "//&
                    "position_copier type unknown.")
        end select
    end subroutine create_position

    subroutine create_orientation(orientation_copier, random_orientation, have_orientations)
        class(Abstract_Coordinates_Copier), allocatable, intent(out) :: orientation_copier
        class(Abstract_Random_Coordinates), intent(in) :: random_orientation
        logical, intent(in) :: have_orientations(:)

        if (any(have_orientations)) then
            allocate(Random_Filling_Coordinates_Copier :: orientation_copier)
        else
            allocate(Null_Coordinates_Copier :: orientation_copier)
        end if
        select type (orientation_copier)
            type is (Random_Filling_Coordinates_Copier)
                call orientation_copier%construct(random_orientation, have_orientations)
            type is (Null_Coordinates_Copier)
            class default
                call error_exit("procedures_coordinates_copier_factory: create_orientation: "//&
                    "orientation_copier type unknown.")
        end select
    end subroutine create_orientation

    subroutine destroy(position_copier)
        class(Abstract_Coordinates_Copier), allocatable, intent(inout) :: position_copier

        if (allocated(position_copier)) then
            call position_copier%destroy()
            deallocate(position_copier)
        end if
    end subroutine destroy

end module procedures_coordinates_copier_factory
