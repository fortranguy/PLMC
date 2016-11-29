module procedures_visitable_walls_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_floor_penetration, only: Abstract_Floor_Penetration
use classes_visitable_walls, only: Abstract_Visitable_Walls, Concrete_Visitable_Walls, &
    Null_Visitable_Walls
use procedures_environment_inquirers, only: use_walls
use classes_min_distance, only: Abstract_Min_Distance

implicit none

private
public :: create, destroy

contains

    subroutine create(visitable_walls, periodic_boxes, floor_penetration, min_distance, &
        generating_data, prefix)
        class(Abstract_Visitable_Walls), allocatable, intent(out) :: visitable_walls(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        class(Abstract_Floor_Penetration), intent(in) :: floor_penetration
        class(Abstract_Min_Distance), intent(in) :: min_distance
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: i_box
        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: gap_centers

        if (use_walls(generating_data, prefix)) then
            allocate(Concrete_Visitable_Walls :: visitable_walls(size(periodic_boxes)))
            data_field = prefix//"Walls.gap between centers"
            call generating_data%get(data_field, gap_centers, data_found)
            call check_data_found(data_field, data_found)
        else
            allocate(Null_Visitable_Walls :: visitable_walls(size(periodic_boxes)))
        end if
        do i_box = 1, size(visitable_walls)
            call visitable_walls(i_box)%construct(periodic_boxes(i_box), gap_centers, &
                floor_penetration, min_distance)
        end do
    end subroutine create

    subroutine destroy(visitable_walls)
        class(Abstract_Visitable_Walls), allocatable, intent(inout) :: visitable_walls(:)

        integer :: i_box

        if (allocated(visitable_walls)) then
            do i_box = size(visitable_walls), 1, -1
                call visitable_walls(i_box)%destroy()
            end do
            deallocate(visitable_walls)
        end if
    end subroutine destroy

end module procedures_visitable_walls_factory
