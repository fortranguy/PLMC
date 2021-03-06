module procedures_density_explorer_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_input_prefixes, only: density_prefix
use json_module, only: json_file
use procedures_errors, only: error_exit, warning_continue
use procedures_checks, only: check_data_found
use classes_number_to_string, only: Concrete_Number_to_String
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_visitable_walls, only: Abstract_Visitable_Walls
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_boxes_factory, only: boxes_create => create, boxes_destroy => destroy
use procedures_environment_inquirers, only: periodicity_is_xy
use classes_density_explorer, only: Abstract_Density_Explorer, Plain_Density_Explorer, &
    Z_Density_Explorer, XZ_Density_Explorer

implicit none

private
public :: create, destroy

contains

    subroutine create(density_explorer, periodic_boxes, max_box_size, visitable_walls, i_component,&
        num_snaps, exploring_data)
        class(Abstract_Density_Explorer), allocatable, intent(out) :: density_explorer
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        real(DP), intent(in) :: max_box_size(:)
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls(:)
        integer, intent(in) :: i_component, num_snaps
        type(json_file) :: exploring_data

        character(len=:), allocatable:: density_type

        call set_density_type(density_type, exploring_data)
        call allocate(density_explorer, density_type)
        call construct(density_explorer, periodic_boxes, max_box_size, visitable_walls, &
            i_component, num_snaps, exploring_data)
    end subroutine create

    subroutine set_density_type(density_type, exploring_data)
        character(len=:), allocatable, intent(out) :: density_type
        type(json_file), intent(inout) :: exploring_data

        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = density_prefix//"name"
        call exploring_data%get(data_field, density_type, data_found)
        call check_data_found(data_field, data_found)
    end subroutine set_density_type

    subroutine allocate(density_explorer, density_type)
        class(Abstract_Density_Explorer), allocatable, intent(out) :: density_explorer
        character(len=*), intent(in):: density_type

        select case(density_type)
            case ("plain")
                allocate(Plain_Density_Explorer :: density_explorer)
            case ("z")
                allocate(Z_Density_Explorer :: density_explorer)
            case ("xz")
                allocate(XZ_Density_Explorer :: density_explorer)
            case default
                call error_exit("procedures_density_explorer_factory: allocate: "//density_type//&
                    " unknown. Please choose between 'plain', 'z' and 'xz'.")
        end select
    end subroutine allocate

    subroutine construct(density_explorer, periodic_boxes, max_box_size, visitable_walls, &
        i_component, num_snaps, exploring_data)
        class(Abstract_Density_Explorer), intent(inout) :: density_explorer
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        real(DP), intent(in) :: max_box_size(:)
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls(:)
        integer, intent(in) :: i_component, num_snaps
        type(json_file) :: exploring_data

        class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domains(:)
        logical :: create_domain
        character(len=:), allocatable ::filename
        type(Concrete_Number_to_String) :: string
        character(len=:), allocatable :: data_field
        logical :: data_found

        create_domain = .true.
        call boxes_create(parallelepiped_domains, periodic_boxes, visitable_walls, create_domain, &
            exploring_data, density_prefix)
        select type (density_explorer)
            type is (Plain_Density_Explorer)
                call density_explorer%construct(parallelepiped_domains(1), num_snaps)
            type is (Z_Density_Explorer)
                call check_periodicity_xy(periodic_boxes)
                block
                    real(DP) :: delta_z
                    data_field = density_prefix//"delta"
                    call exploring_data%get(data_field, delta_z, data_found)
                    call check_data_found(data_field, data_found)
                    filename = "z_density_"//string%get(i_component)//".out"
                    call density_explorer%construct(max_box_size, parallelepiped_domains(1), &
                        delta_z, num_snaps, filename)
                end block
            type is (XZ_Density_Explorer)
                call check_periodicity_xy(periodic_boxes)
                block
                    real(DP) :: delta_xz
                    data_field = density_prefix//"delta"
                    call exploring_data%get(data_field, delta_xz, data_found)
                    call check_data_found(data_field, data_found)
                    filename = "xz_density_"//string%get(i_component)//".out"
                    call density_explorer%construct(max_box_size, parallelepiped_domains(1), &
                        delta_xz, num_snaps, filename)
                end block
            class default
                call error_exit("procedures_density_explorer_factory: construct: "//&
                    "density_explorer type unknown.")
        end select
        call boxes_destroy(parallelepiped_domains)
    end subroutine construct

    subroutine check_periodicity_xy(periodic_boxes)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)

        if (all(.not.periodicity_is_xy(periodic_boxes))) then
            call warning_continue("Periodicity is not XY.")
        end if
    end subroutine check_periodicity_xy

    subroutine destroy(density_explorer)
        class(Abstract_Density_Explorer), allocatable, intent(inout) :: density_explorer

        if (allocated(density_explorer)) then
            call density_explorer%destroy()
            deallocate(density_explorer)
        end if
    end subroutine destroy

end module procedures_density_explorer_factory
