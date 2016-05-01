module procedures_visitable_walls_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_floor_penetration, only: Abstract_Floor_Penetration
use classes_visitable_walls, only: Abstract_Visitable_Walls, Concrete_Visitable_Walls, &
    Null_Visitable_Walls
use procedures_property_inquirers, only: use_walls

implicit none

private
public :: create, destroy

contains

    subroutine create(walls, periodic_box, floor_penetration, input_data, prefix)
        class(Abstract_Visitable_Walls), allocatable, intent(out) :: walls
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Floor_Penetration), intent(in) :: floor_penetration
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: gap

        if (use_walls(input_data, prefix)) then
            allocate(Concrete_Visitable_Walls :: walls)
            data_field = prefix//"Walls.gap"
            call input_data%get(data_field, gap, data_found)
            call check_data_found(data_field, data_found)
        else
            allocate(Null_Visitable_Walls :: walls)
        end if
        call walls%construct(periodic_box, gap, floor_penetration)
    end subroutine create

    subroutine destroy(walls)
        class(Abstract_Visitable_Walls), allocatable, intent(inout) :: walls

        if (allocated(walls)) then
            call walls%destroy()
            deallocate(walls)
        end if
    end subroutine destroy

end module procedures_visitable_walls_factory