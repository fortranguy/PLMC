module procedures_random_coordinates_factory

use json_module, only: json_file
use procedures_errors, only: error_exit
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_random_coordinates, only: Abstract_Random_Coordinates, Null_Random_Coordinates
use classes_random_position, only: Concrete_Random_Position
use classes_random_orientation, only: Concrete_Random_Orientation

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_random_positions
    module procedure :: create_random_orientation
end interface create

interface destroy
    module procedure :: destroy_line
    module procedure :: destroy_element
end interface destroy

contains

    subroutine create_random_positions(random_positions, parallelepiped_domains, have_positions, &
        can_exchange)
        class(Abstract_Random_Coordinates), allocatable, intent(out) :: random_positions(:)
        class(Abstract_Parallelepiped_Domain), intent(in) :: parallelepiped_domains(:)
        logical, intent(in) :: have_positions(:, :)
        logical, intent(in) :: can_exchange(:)

        integer :: i_box

        if (any(have_positions) .and. any(can_exchange)) then
            allocate(Concrete_Random_Position :: random_positions(size(parallelepiped_domains)))
        else
            allocate(Null_Random_Coordinates :: random_positions(size(parallelepiped_domains)))
        end if
        select type(random_positions)
            type is (Concrete_Random_Position)
                do i_box = 1, size(random_positions)
                    call random_positions(i_box)%construct(parallelepiped_domains(i_box), &
                        have_positions(:, i_box))
                end do
            type is (Null_Random_Coordinates)
            class default
                call error_exit("procedures_random_coordinates_factory: create_random_positions: "&
                    "random_positions: unknown type.")
        end select
    end subroutine create_random_positions

    subroutine create_random_orientation(random_orientation, have_orientations, can_exchange)
        class(Abstract_Random_Coordinates), allocatable, intent(out) :: random_orientation
        logical, dimension(:), intent(in) :: have_orientations, can_exchange

        if (any(have_orientations) .and. any(can_exchange)) then
            allocate(Concrete_Random_Orientation :: random_orientation)
        else
            allocate(Null_Random_Coordinates :: random_orientation)
        end if
        select type(random_orientation)
            type is (Concrete_Random_Orientation)
                call random_orientation%construct(have_orientations)
            type is (Null_Random_Coordinates)
            class default
                call error_exit("procedures_random_coordinates_factory: create_random_orientation:"&
                    " random_orientation: unknown type.")
        end select
    end subroutine create_random_orientation

    subroutine destroy_line(random_coordinates)
        class(Abstract_Random_Coordinates), allocatable, intent(inout) :: random_coordinates(:)

        integer :: i_box

        if (allocated(random_coordinates)) then
            do i_box = size(random_coordinates), 1, -1
                call random_coordinates(i_box)%destroy()
            end do
            deallocate(random_coordinates)
        end if
    end subroutine destroy_line

    subroutine destroy_element(random_coordinates)
        class(Abstract_Random_Coordinates), allocatable, intent(inout) :: random_coordinates

        if (allocated(random_coordinates)) then
            call random_coordinates%destroy()
            deallocate(random_coordinates)
        end if
    end subroutine destroy_element

end module procedures_random_coordinates_factory
