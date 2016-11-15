program density

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use data_input_prefixes, only: environment_prefix, mixture_prefix, density_prefix
use data_arguments, only: i_generating, i_exploring, num_json_arguments
use json_module, only: json_file
use procedures_json_data_factory, only: json_data_create => create, json_data_destroy => destroy
use procedures_checks, only: check_in_range
use procedures_command_arguments, only: create_filename_from_argument
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_boxes_factory, only: boxes_create => create, boxes_destroy => destroy
use classes_floor_penetration, only: Abstract_Floor_Penetration
use classes_visitable_walls, only: Abstract_Visitable_Walls
use procedures_walls_factory, only: walls_create => create, walls_destroy => destroy
use procedures_environment_inquirers, only: use_walls
use classes_min_distance, only: Abstract_Min_Distance
use procedures_hard_core_factory, only: hard_core_create => create, hard_core_destroy => destroy
use types_raw_coordinates, only: Concrete_Raw_Coordinates
use procedures_mixture_inquirers, only: property_num_components => num_components, &
    property_i_component => i_component
use types_component_coordinates_reader_selector, only: Component_Coordinates_Reader_Selector
use procedures_complete_coordinates_reader, only: complete_coordinates_read => read, &
    complete_coordinates_deallocate => deallocate
use procedures_plmc_iterations, only: plmc_set_num_snaps
use classes_density_explorer, only: Abstract_Density_Explorer
use procedures_density_explorer_factory, only: density_explorer_create => create, &
    density_explorer_destroy => destroy
use procedures_plmc_help, only: plmc_catch_density_help

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_boxes(:)
    real(DP), dimension(num_dimensions) :: max_box_size, box_size
    class(Abstract_Visitable_Walls), allocatable :: visitable_walls(:)
    class(Abstract_Floor_Penetration), allocatable :: floor_penetration
    class(Abstract_Min_Distance), allocatable :: wall_min_distance
    integer :: num_components, i_component
    type(Component_Coordinates_Reader_Selector) :: selector
    integer :: num_snaps, i_snap
    class(Abstract_Density_Explorer), allocatable :: density_explorer
    character(len=:), allocatable :: snap_filename
    type(Concrete_Raw_Coordinates) :: raw_coordinates
    type(json_file) :: generating_data, exploring_data

    call plmc_catch_density_help()
    call json_data_create(generating_data, i_generating)
    call plmc_set_num_snaps(num_snaps, generating_data)
    call boxes_create(periodic_boxes, generating_data, environment_prefix, unique=.true.)
    call walls_create(floor_penetration, generating_data, environment_prefix)
    call hard_core_create(wall_min_distance, use_walls(floor_penetration), generating_data, &
        environment_prefix//"Walls.")
    call walls_create(visitable_walls, periodic_boxes, floor_penetration, wall_min_distance, &
        generating_data, environment_prefix)
    num_components = property_num_components(generating_data, mixture_prefix)
    call json_data_destroy(generating_data)

    max_box_size = 0._DP
    do i_snap = 1, num_snaps
        call create_filename_from_argument(snap_filename, num_json_arguments + i_snap)
        call complete_coordinates_read(box_size, snap_filename)
        if (any(max_box_size < box_size)) max_box_size = box_size
        if (allocated(snap_filename)) deallocate(snap_filename)
    end do

    call json_data_create(exploring_data, i_exploring)
    i_component = property_i_component(exploring_data, density_prefix)
    call check_in_range("density", num_components, "i_component", i_component)
    call density_explorer_create(density_explorer, periodic_boxes, max_box_size, visitable_walls, &
        i_component, num_snaps, exploring_data)
    call json_data_destroy(exploring_data)

    selector%read_positions = .true.
    selector%read_orientations = .false.
    do i_snap = 1, num_snaps
        call create_filename_from_argument(snap_filename, num_json_arguments + i_snap)
        call complete_coordinates_read(box_size, raw_coordinates, num_components, i_component, &
            snap_filename, selector)
        call periodic_boxes(1)%set(box_size)
        call density_explorer%fill(i_snap, raw_coordinates%positions)
        call complete_coordinates_deallocate(raw_coordinates)
        if (allocated(snap_filename)) deallocate(snap_filename)
    end do
    call density_explorer%write()

    call density_explorer_destroy(density_explorer)
    call walls_destroy(visitable_walls)
    call hard_core_destroy(wall_min_distance)
    call walls_destroy(floor_penetration)
    call boxes_destroy(periodic_boxes)

end program density
