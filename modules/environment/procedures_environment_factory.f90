module procedures_environment_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: warning_continue
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_box_factory, only: box_create => create, box_destroy => destroy
use procedures_temperature_factory, only: temperature_create => create, &
    temperature_destroy => destroy
use classes_field_expression, only: Abstract_Field_Expression
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_field_factory, only: field_create => create, field_destroy => destroy
use procedures_permittivity_factory, only: permittivity_create => create, &
    permittivity_destroy => destroy
use classes_floor_penetration, only: Abstract_Floor_Penetration
use classes_walls_potential, only: Abstract_Walls_Potential
use procedures_walls_factory, only: walls_create => create, walls_destroy => destroy
use procedures_component_factory, only: component_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use procedures_property_inquirers, only: periodicity_is_xyz, periodicity_is_xy, &
    apply_external_field, use_walls

implicit none

private
public :: environment_create, environment_destroy

contains

    subroutine environment_create(environment, input_data, prefix)
        type(Environment_Wrapper), intent(out) :: environment
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        class(Abstract_Field_Expression), allocatable :: field_expression
        class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
        class(Abstract_Floor_Penetration), allocatable :: floor_penetration
        logical :: field_applied

        call box_create(environment%periodic_box, input_data, prefix)
        call temperature_create(environment%temperature, input_data, prefix)
        field_applied = apply_external_field(input_data, prefix)
        call permittivity_create(environment%permittivity, input_data, prefix)
        call field_create(field_expression, environment%permittivity, field_applied, &
            input_data, prefix)
        call box_create(parallelepiped_domain, environment%periodic_box, field_applied, &
            input_data, prefix//"External Field.")
        call field_create(environment%external_field, parallelepiped_domain, field_expression, &
            field_applied)
        call field_destroy(field_expression)
        call box_destroy(parallelepiped_domain)
        call box_create(environment%reciprocal_lattice, environment%periodic_box, &
            input_data, prefix)
        call walls_create(floor_penetration, input_data, prefix)
        call walls_create(environment%walls_potential, environment%periodic_box, &
            floor_penetration, input_data, prefix)
        call walls_destroy(floor_penetration)
        call box_create(environment%box_size_checker, environment%reciprocal_lattice, &
            environment%walls_potential)

        call environment_check(environment%periodic_box, environment%walls_potential)
        call environment%box_size_checker%check()
    end subroutine environment_create

    subroutine environment_destroy(environment)
        type(Environment_Wrapper), intent(inout) :: environment

        call box_destroy(environment%box_size_checker)
        call walls_destroy(environment%walls_potential)
        call box_destroy(environment%reciprocal_lattice)
        call field_destroy(environment%external_field)
        call permittivity_destroy(environment%permittivity)
        call temperature_destroy(environment%temperature)
        call box_destroy(environment%periodic_box)
    end subroutine environment_destroy

    subroutine environment_check(periodic_box, walls_potential)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Walls_Potential), intent(in) :: walls_potential

        if (periodicity_is_xyz(periodic_box) .and. use_walls(walls_potential)) then
            call warning_continue("environment_check: periodicity is XYZ but walls are used.")
        end if
        if (periodicity_is_xy(periodic_box) .and. .not.use_walls(walls_potential)) then
            call warning_continue("environment_check: periodicity is XY but walls are not used.")
        end if
    end subroutine environment_check

end module procedures_environment_factory
