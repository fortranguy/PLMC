module procedures_environment_inquirers

use json_module, only: json_file
use procedures_checks, only: check_data_found
use procedures_property_inquirers, only: logical_from_json
use classes_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain, Null_Parallelepiped_Domain
use classes_permittivity, only: Abstract_Permittivity, Concrete_Permittivity
use classes_beta_pressure, only: Abstract_Beta_Pressure, Concrete_Beta_Pressure
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice, Concrete_Reciprocal_Lattice
use classes_external_field, only: Abstract_External_Field, Null_External_Field
use classes_floor_penetration, only: Abstract_Floor_Penetration, Null_Floor_Penetration
use classes_visitable_walls, only: Abstract_Visitable_Walls, Concrete_Visitable_Walls
use classes_changed_box_size, only: Abstract_Changed_Box_Size, Concrete_Changed_Box_Size

implicit none

private
public :: num_boxes, periodicity_is_xyz, periodicity_is_xy, total_volume_can_change, &
    box_size_can_change, use_domain, apply_external_field, use_permittivity, &
    use_reciprocal_lattice, use_walls

interface num_boxes
    module procedure :: num_boxes_from_json
end interface num_boxes

interface total_volume_can_change
    module procedure :: total_volume_can_change_from_json
    module procedure :: total_volume_can_change_from_beta_pressure
end interface total_volume_can_change

interface box_size_can_change
    module procedure :: box_size_can_change_from_changed
end interface box_size_can_change

interface apply_external_field
    module procedure :: apply_external_field_from_json
    module procedure :: apply_external_field_from
end interface apply_external_field

interface use_permittivity
    module procedure :: use_permittivity_from_json
    module procedure :: use_permittivity_from
end interface

interface use_reciprocal_lattice
    module procedure :: use_reciprocal_lattice_from_json
    module procedure :: use_reciprocal_lattice_from
end interface use_reciprocal_lattice

interface use_walls
    module procedure :: use_walls_from_json
    module procedure :: use_walls_from_floor_penetration
    module procedure :: use_walls_from_walls
end interface use_walls

contains

    integer function num_boxes_from_json(generating_data, prefix) result(num_boxes)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = prefix//"Boxes.number"
        call generating_data%get(data_field, num_boxes, data_found)
        call check_data_found(data_field, data_found)
    end function num_boxes_from_json

    elemental logical function periodicity_is_xyz(periodic_box)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        select type (periodic_box)
            type is (XYZ_Periodic_Box)
                periodicity_is_xyz = .true.
            class default
                periodicity_is_xyz = .false.
        end select
    end function periodicity_is_xyz

    elemental logical function periodicity_is_xy(periodic_box)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        select type (periodic_box)
            type is (XY_Periodic_Box)
                periodicity_is_xy = .true.
            class default
                periodicity_is_xy = .false.
        end select
    end function periodicity_is_xy

    logical function total_volume_can_change_from_json(generating_data, prefix) &
        result(total_volume_can_change)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        total_volume_can_change = logical_from_json(generating_data, prefix//&
            "Boxes.total volume can change")
    end function total_volume_can_change_from_json

    pure logical function total_volume_can_change_from_beta_pressure(beta_pressure) &
        result(total_volume_can_change)
        class(Abstract_Beta_Pressure), intent(in) :: beta_pressure

        select type (beta_pressure)
            type is (Concrete_Beta_Pressure)
                total_volume_can_change = .true.
            class default
                total_volume_can_change = .false.
        end select
    end function total_volume_can_change_from_beta_pressure

    elemental logical function box_size_can_change_from_changed(changed_box_size) &
        result(box_size_can_change)
        class(Abstract_Changed_Box_Size), intent(in) :: changed_box_size

        select type (changed_box_size)
            type is (Concrete_Changed_Box_Size)
                box_size_can_change = .true.
            class default
                box_size_can_change = .false.
        end select
    end function box_size_can_change_from_changed

    elemental logical function use_domain(parallelepiped_domain)
        class(Abstract_Parallelepiped_Domain), intent(in) :: parallelepiped_domain

        select type (parallelepiped_domain)
            type is (Null_Parallelepiped_Domain)
                use_domain = .false.
            class default
                use_domain = .true.
        end select
    end function use_domain

    logical function apply_external_field_from_json(generating_data, prefix) &
        result(apply_external_field)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        apply_external_field = logical_from_json(generating_data, prefix//"External Field.apply")
    end function apply_external_field_from_json

    pure logical function apply_external_field_from(external_field) result(apply_external_field)
        class(Abstract_External_Field), intent(in) :: external_field

        select type (external_field)
            type is (Null_External_Field)
                apply_external_field = .false.
            class default
                apply_external_field = .true.
        end select
    end function apply_external_field_from

    logical function use_permittivity_from_json(generating_data, prefix) result(use_permittivity)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        use_permittivity = logical_from_json(generating_data, prefix//"Permittivity.use")
    end function use_permittivity_from_json

    pure logical function use_permittivity_from(permittivity) result(use_permittivity)
        class(Abstract_Permittivity), intent(in) :: permittivity

        select type (permittivity)
            type is (Concrete_Permittivity)
                use_permittivity = .true.
            class default
                use_permittivity = .false.
        end select
    end function use_permittivity_from

    logical function use_reciprocal_lattice_from_json(generating_data, prefix) &
        result(use_reciprocal_lattice)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        use_reciprocal_lattice = logical_from_json(generating_data, &
            prefix//"Reciprocal Lattice.use")
    end function use_reciprocal_lattice_from_json

    elemental logical function use_reciprocal_lattice_from(reciprocal_lattice) &
        result(use_reciprocal_lattice)
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice

        select type (reciprocal_lattice)
            type is (Concrete_Reciprocal_Lattice)
                use_reciprocal_lattice = .true.
            class default
                use_reciprocal_lattice = .false.
        end select
    end function use_reciprocal_lattice_from

    logical function use_walls_from_json(generating_data, prefix) result(use_walls)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        use_walls = logical_from_json(generating_data, prefix//"Walls.use")
    end function use_walls_from_json

    pure logical function use_walls_from_floor_penetration(floor_penetration) result(use_walls)
        class(Abstract_Floor_Penetration), intent(in) :: floor_penetration

        select type (floor_penetration)
            type is (Null_Floor_Penetration)
                use_walls = .false.
            class default
                use_walls = .true.
        end select
    end function use_walls_from_floor_penetration

    elemental logical function use_walls_from_walls(visitable_walls) result(use_walls)
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls

        select type (visitable_walls)
            type is (Concrete_Visitable_Walls)
                use_walls = .true.
            class default
                use_walls = .false.
        end select
    end function use_walls_from_walls

end module procedures_environment_inquirers
