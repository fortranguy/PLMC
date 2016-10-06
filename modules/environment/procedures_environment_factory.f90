module procedures_environment_factory

use json_module, only: json_file
use procedures_errors, only: error_exit, warning_continue
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_size_checker, only: Abstract_Box_Size_Checker
use procedures_boxes_factory, only: boxes_create => create, boxes_destroy => destroy
use procedures_beta_pressure_factory, only: beta_pressure_create => create, &
    beta_pressure_destroy => destroy
use procedures_temperature_factory, only: temperature_create => create, &
    temperature_destroy => destroy
use classes_field_expression, only: Abstract_Field_Expression
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_fields_factory, only: fields_create => create, fields_destroy => destroy
use procedures_permittivity_factory, only: permittivity_create => create, &
    permittivity_destroy => destroy
use classes_floor_penetration, only: Abstract_Floor_Penetration
use classes_visitable_walls, only: Abstract_Visitable_Walls
use procedures_walls_factory, only: walls_create => create, walls_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_inquirers, only: periodicity_is_xyz, periodicity_is_xy, &
    total_volume_can_change, apply_external_field, use_walls
use procedures_hard_core_factory, only: hard_core_create => create, hard_core_destroy => destroy

implicit none

private
public :: environment_create, environment_destroy

contains

    subroutine environment_create(environment, generating_data, prefix)
        type(Environment_Wrapper), intent(out) :: environment
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        class(Abstract_Field_Expression), allocatable :: field_expression
        class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domains(:)
        class(Abstract_Floor_Penetration), allocatable :: floor_penetration
        logical :: field_applied

        call boxes_create(environment%periodic_boxes, generating_data, prefix)
        call beta_pressure_create(environment%beta_pressure, &
            total_volume_can_change(generating_data, prefix), generating_data, prefix)
        call temperature_create(environment%temperature, generating_data, prefix)
        field_applied = apply_external_field(generating_data, prefix)
        call permittivity_create(environment%permittivity, generating_data, prefix)
        call walls_create(floor_penetration, generating_data, prefix)
        call hard_core_create(environment%wall_min_distance, use_walls(floor_penetration), &
            generating_data, prefix//"Walls.")
        call walls_create(environment%gemc_visitable_walls, environment%periodic_boxes, &
            floor_penetration,environment%wall_min_distance, generating_data, prefix)
        call walls_destroy(floor_penetration)
        call fields_create(field_expression, environment%permittivity, field_applied, &
            generating_data, prefix)
        call boxes_create(parallelepiped_domains, environment%periodic_boxes, environment%&
            gemc_visitable_walls, field_applied, generating_data, prefix//"External Field.")
        call fields_create(environment%external_fields, parallelepiped_domains, field_expression, &
            field_applied)
        call fields_destroy(field_expression)
        call boxes_destroy(parallelepiped_domains)
        call boxes_create(environment%reciprocal_lattices, environment%periodic_boxes, &
            generating_data, prefix)
        call boxes_create(environment%box_size_checkers, environment%reciprocal_lattices, &
            environment%gemc_visitable_walls)
        if (all(periodicity_is_xyz(environment%periodic_boxes))) then
            call boxes_create(environment%accessible_domains, environment%periodic_boxes, &
                needed=.true.)
        else if (all(periodicity_is_xy(environment%periodic_boxes))) then
            call boxes_create(environment%accessible_domains, environment%periodic_boxes, &
                environment%gemc_visitable_walls, needed=.true.)
        else
            call error_exit("procedures_environment_factory: environment_create: "//&
                "box periodicity is unknown.")
        end if

        call check(environment%periodic_boxes, environment%gemc_visitable_walls, environment%&
            box_size_checkers)
    end subroutine environment_create

    subroutine check(periodic_boxes, visitable_walls, box_size_checkers)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls(:)
        class(Abstract_Box_Size_Checker), intent(in) :: box_size_checkers(:)

        integer :: i_box

        if (all(periodicity_is_xyz(periodic_boxes)) .and. all(use_walls(visitable_walls))) then
            call warning_continue("procedures_environment_factory: check: "//&
                "periodicity is XYZ but walls are used.")
        end if
        if (all(periodicity_is_xy(periodic_boxes)) .and. all(.not.use_walls(visitable_walls))) then
            call warning_continue("procedures_environment_factory: check: "//&
                "periodicity is XY but walls are not used.")
        end if

        do i_box = 1, size(box_size_checkers)
            call box_size_checkers(i_box)%check()
        end do
    end subroutine check

    subroutine environment_destroy(environment)
        type(Environment_Wrapper), intent(inout) :: environment

        call boxes_destroy(environment%accessible_domains)
        call boxes_destroy(environment%box_size_checkers)
        call boxes_destroy(environment%reciprocal_lattices)
        call fields_destroy(environment%external_fields)
        call walls_destroy(environment%gemc_visitable_walls)
        call hard_core_destroy(environment%wall_min_distance)
        call permittivity_destroy(environment%permittivity)
        call temperature_destroy(environment%temperature)
        call beta_pressure_destroy(environment%beta_pressure)
        call boxes_destroy(environment%periodic_boxes)
    end subroutine environment_destroy

end module procedures_environment_factory
