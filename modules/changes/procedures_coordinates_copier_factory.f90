module procedures_coordinates_copier_factory

use procedures_errors, only: error_exit
use classes_random_coordinates, only: Abstract_Random_Coordinates
use classes_coordinates_copier, only: Abstract_Coordinates_Copier, &
    Random_Filling_Coordinates_Copier, Null_Coordinates_Copier

implicit none

private
public :: create_position, create_orientation, destroy

interface destroy
    module procedure :: destroy_line
    module procedure :: destroy_element
end interface destroy

contains

    subroutine create_position(position_copiers, random_positions, have_positions, can_exchange)
        class(Abstract_Coordinates_Copier), allocatable, intent(out) :: position_copiers(:)
        class(Abstract_Random_Coordinates), intent(in) :: random_positions(:)
        logical, intent(in) :: have_positions(:, :)
        logical, intent(in) :: can_exchange(:)

        integer :: i_box

        if (any(have_positions) .and. any(can_exchange)) then
            allocate(Random_Filling_Coordinates_Copier :: position_copiers(size(random_positions)))
        else
            allocate(Null_Coordinates_Copier :: position_copiers(size(random_positions)))
        end if
        select type (position_copiers)
            type is (Random_Filling_Coordinates_Copier)
                do i_box = 1, size(position_copiers)
                    call position_copiers(i_box)%construct(random_positions(i_box), &
                        have_positions(:, i_box))
                end do
            type is (Null_Coordinates_Copier)
            class default
                call error_exit("procedures_coordinates_copier_factory: create_position: "//&
                    "position_copiers type unknown.")
        end select
    end subroutine create_position

    subroutine create_orientation(orientation_copier, random_orientation, have_orientations, &
        can_exchange)
        class(Abstract_Coordinates_Copier), allocatable, intent(out) :: orientation_copier
        class(Abstract_Random_Coordinates), intent(in) :: random_orientation
        logical, intent(in) :: have_orientations(:)
        logical, intent(in) :: can_exchange(:)

        if (any(have_orientations) .and. any(can_exchange)) then
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

    subroutine destroy_line(coordinates_copiers)
        class(Abstract_Coordinates_Copier), allocatable, intent(inout) :: coordinates_copiers(:)

        integer :: i_box

        if (allocated(coordinates_copiers)) then
            do i_box = size(coordinates_copiers), 1, -1
                call coordinates_copiers(i_box)%destroy()
            end do
            deallocate(coordinates_copiers)
        end if
    end subroutine destroy_line

    subroutine destroy_element(coordinates_copier)
        class(Abstract_Coordinates_Copier), allocatable, intent(inout) :: coordinates_copier

        if (allocated(coordinates_copier)) then
            call coordinates_copier%destroy()
            deallocate(coordinates_copier)
        end if
    end subroutine destroy_element

end module procedures_coordinates_copier_factory
