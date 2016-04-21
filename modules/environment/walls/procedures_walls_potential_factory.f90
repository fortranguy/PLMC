module procedures_walls_potential_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_floor_penetration, only: Abstract_Floor_Penetration
use classes_walls_potential, only: Abstract_Walls_Potential, Concrete_Walls_Potential, &
    Null_Walls_Potential
use procedures_property_inquirers, only: use_walls

implicit none

private
public :: create, destroy

contains

    subroutine create(walls_potential, periodic_box, floor_penetration, input_data, prefix)
        class(Abstract_Walls_Potential), allocatable, intent(out) :: walls_potential
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Floor_Penetration), intent(in) :: floor_penetration
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: gap

        if (use_walls(input_data, prefix)) then
            allocate(Concrete_Walls_Potential :: walls_potential)
            data_field = prefix//"Walls.gap"
            call input_data%get(data_field, gap, data_found)
            call check_data_found(data_field, data_found)
        else
            allocate(Null_Walls_Potential :: walls_potential)
        end if
        call walls_potential%construct(periodic_box, gap, floor_penetration)
    end subroutine create

    subroutine destroy(walls_potential)
        class(Abstract_Walls_Potential), allocatable, intent(inout) :: walls_potential

        if (allocated(walls_potential)) then
            call walls_potential%destroy()
            deallocate(walls_potential)
        end if
    end subroutine destroy

end module procedures_walls_potential_factory
