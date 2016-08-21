module procedures_environment_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit, warning_continue
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_box_factory, only: box_create => create, box_destroy => destroy
use procedures_beta_pressure_factory, only: beta_pressure_create => create, &
    beta_pressure_destroy => destroy
use procedures_temperature_factory, only: temperature_create => create, &
    temperature_destroy => destroy
use classes_field_expression, only: Abstract_Field_Expression
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_field_factory, only: field_create => create, field_destroy => destroy
use procedures_permittivity_factory, only: permittivity_create => create, &
    permittivity_destroy => destroy
use classes_floor_penetration, only: Abstract_Floor_Penetration
use classes_visitable_walls, only: Abstract_Visitable_Walls
use procedures_walls_factory, only: walls_create => create, walls_destroy => destroy
use procedures_component_factory, only: component_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use procedures_hard_core_factory, only: hard_core_create => create, hard_core_destroy => destroy
use procedures_property_inquirers, only: periodicity_is_xyz, periodicity_is_xy, &
    box_size_can_change, apply_external_field, use_walls

implicit none

private
public :: environment_create, environment_destroy

contains

    subroutine environment_create(environment, generating_data, prefix)
        type(Environment_Wrapper), intent(out) :: environment
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        class(Abstract_Field_Expression), allocatable :: field_expression
        class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
        class(Abstract_Floor_Penetration), allocatable :: floor_penetration
        logical :: field_applied

        call box_create(environment%periodic_box, generating_data, prefix)
        call beta_pressure_create(environment%beta_pressure, box_size_can_change(generating_data, &
            prefix), generating_data, prefix)
        call temperature_create(environment%temperature, generating_data, prefix)
        field_applied = apply_external_field(generating_data, prefix)
        call permittivity_create(environment%permittivity, generating_data, prefix)
        call walls_create(floor_penetration, generating_data, prefix)
        call hard_core_create(environment%wall_min_distance, use_walls(floor_penetration), &
            generating_data, prefix//"Walls.")
        call walls_create(environment%visitable_walls, environment%periodic_box, floor_penetration,&
            environment%wall_min_distance, generating_data, prefix)
        call walls_destroy(floor_penetration)
        call field_create(field_expression, environment%permittivity, field_applied, &
            generating_data, prefix)
        call box_create(parallelepiped_domain, environment%periodic_box, environment%&
            visitable_walls, field_applied, generating_data, prefix//"External Field.")
        call field_create(environment%external_field, parallelepiped_domain, field_expression, &
            field_applied)
        call field_destroy(field_expression)
        call box_destroy(parallelepiped_domain)
        call box_create(environment%reciprocal_lattice, environment%periodic_box, &
            generating_data, prefix)
        call box_create(environment%box_size_checker, environment%reciprocal_lattice, &
            environment%visitable_walls)
        if (periodicity_is_xyz(environment%periodic_box)) then
            call box_create(environment%accessible_domain, environment%periodic_box, .true.)
        else if (periodicity_is_xy(environment%periodic_box)) then
            call box_create(environment%accessible_domain, environment%periodic_box, environment%&
                visitable_walls, .true.)
        else
            call error_exit("procedures_environment_factory: environment_create: "//&
                "box periodicity is unknown.")
        end if

        call check(environment%periodic_box, environment%visitable_walls)
        call environment%box_size_checker%check()
    end subroutine environment_create

    subroutine environment_destroy(environment)
        type(Environment_Wrapper), intent(inout) :: environment

        call box_destroy(environment%accessible_domain)
        call box_destroy(environment%box_size_checker)
        call box_destroy(environment%reciprocal_lattice)
        call field_destroy(environment%external_field)
        call walls_destroy(environment%visitable_walls)
        call hard_core_destroy(environment%wall_min_distance)
        call permittivity_destroy(environment%permittivity)
        call temperature_destroy(environment%temperature)
        call beta_pressure_destroy(environment%beta_pressure)
        call box_destroy(environment%periodic_box)
    end subroutine environment_destroy

    subroutine check(periodic_box, walls)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Visitable_Walls), intent(in) :: walls

        if (periodicity_is_xyz(periodic_box) .and. use_walls(walls)) then
            call warning_continue("procedures_environment_factory: check: "//&
                "periodicity is XYZ but walls are used.")
        end if
        if (periodicity_is_xy(periodic_box) .and. .not.use_walls(walls)) then
            call warning_continue("procedures_environment_factory: check: "//&
                "periodicity is XY but walls are not used.")
        end if
    end subroutine check

end module procedures_environment_factory
