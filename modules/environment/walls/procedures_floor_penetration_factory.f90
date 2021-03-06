module procedures_floor_penetration_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_floor_penetration, only: Abstract_Floor_Penetration, Flat_Floor_Penetration, &
    Centered_Block_Penetration, Null_Floor_Penetration
use procedures_environment_inquirers, only: use_walls

implicit none

private
public :: create, destroy

contains

    subroutine create(floor_penetration, generating_data, prefix)
        class(Abstract_Floor_Penetration), allocatable, intent(out) :: floor_penetration
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        call allocate_floor_penetration(floor_penetration, generating_data, prefix)
        call set_floor_penetration(floor_penetration, generating_data, prefix)
    end subroutine create

    subroutine allocate_floor_penetration(floor_penetration, generating_data, prefix)
        class(Abstract_Floor_Penetration), allocatable, intent(out) :: floor_penetration
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: walls_name, data_field
        logical :: data_found

        if (use_walls(generating_data, prefix)) then
            data_field = prefix//"Walls.name"
            call generating_data%get(data_field, walls_name, data_found)
            call check_data_found(data_field, data_found)
            select case (walls_name)
                case ("flat")
                    allocate(Flat_Floor_Penetration :: floor_penetration)
                case ("blocks")
                    allocate(Centered_Block_Penetration :: floor_penetration)
                case default
                    call error_exit(walls_name//" walls_name unknown. "//&
                        "Choose between 'flat' and 'blocks'.")
            end select
        else
            allocate(Null_Floor_Penetration :: floor_penetration)
        end if
    end subroutine allocate_floor_penetration

    subroutine set_floor_penetration(floor_penetration, generating_data, prefix)
        class(Abstract_Floor_Penetration), intent(inout) :: floor_penetration
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        select type (floor_penetration)
            type is (Flat_Floor_Penetration)
            type is (Centered_Block_Penetration)
                block
                    real(DP), allocatable :: block_size(:)
                    real(DP) ::block_radius
                    data_field = prefix//"Walls.size"
                    call generating_data%get(data_field, block_size, data_found)
                    call check_data_found(data_field, data_found)
                    data_field = prefix//"Walls.radius"
                    call generating_data%get(data_field, block_radius, data_found)
                    call check_data_found(data_field, data_found)
                    call floor_penetration%set(block_size, block_radius)
                end block
            type is (Null_Floor_Penetration)
            class default
                call error_exit("floor_penetration type unknown.")
        end select
    end subroutine set_floor_penetration

    subroutine destroy(floor_penetration)
        class(Abstract_Floor_Penetration), allocatable, intent(inout) :: floor_penetration

        if (allocated(floor_penetration)) deallocate(floor_penetration)
    end subroutine destroy

end module procedures_floor_penetration_factory
