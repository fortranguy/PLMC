module procedures_random_coordinates_factory

use json_module, only: json_file
use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_box_factory, only: box_create => create, box_destroy => destroy
use classes_random_coordinates, only: Abstract_Random_Coordinates, Null_Random_Coordinates
use classes_random_position, only: Concrete_Random_Position
use classes_random_orientation, only: Concrete_Random_Orientation

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_random_position
    module procedure :: create_random_orientation
end interface create

contains

    subroutine create_random_position(random_position, periodic_box, can_exchange, have_positions, &
        json_data, prefix)
        class(Abstract_Random_Coordinates), allocatable, intent(out) :: random_position
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, dimension(:), intent(in) :: can_exchange, have_positions
        type(json_file), intent(inout) :: json_data
        character(len=*), intent(in) :: prefix

        class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain

        if (any(can_exchange) .and. any(have_positions)) then
            allocate(Concrete_Random_Position :: random_position)
        else
            allocate(Null_Random_Coordinates :: random_position)
        end if
        select type(random_position)
            type is (Concrete_Random_Position)
                call box_create(parallelepiped_domain, periodic_box, any(can_exchange) .and. &
                    any(have_positions), json_data, prefix)
                call random_position%construct(periodic_box, parallelepiped_domain, have_positions)
                call box_destroy(parallelepiped_domain)
            type is (Null_Random_Coordinates)
            class default
                call error_exit("procedures_random_coordinates_factory: create_random_position: "&
                    "random_position: unknown type.")
        end select
    end subroutine create_random_position

    subroutine create_random_orientation(random_orientation, can_exchange, have_orientations)
        class(Abstract_Random_Coordinates), allocatable, intent(out) :: random_orientation
        logical, dimension(:), intent(in) :: can_exchange, have_orientations

        if (any(can_exchange) .and. any(have_orientations)) then
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

    subroutine destroy(random_coordinates)
        class(Abstract_Random_Coordinates), allocatable, intent(inout) :: random_coordinates

        if (allocated(random_coordinates)) then
            call random_coordinates%destroy()
            deallocate(random_coordinates)
        end if
    end subroutine destroy

end module procedures_random_coordinates_factory
