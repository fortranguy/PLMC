program density

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_input_prefixes, only: environment_prefix, mixture_prefix, density_prefix
use data_constants, only: num_dimensions
use json_module, only: json_file
use procedures_json_data_factory, only: json_data_create_input => create_input, &
    json_data_destroy_input => destroy_input
use procedures_errors, only: error_exit
use procedures_checks, only: check_in_range
use procedures_command_arguments, only: create_filename_from_argument
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_box_factory, only: box_create => create, box_destroy => destroy
use classes_floor_penetration, only: Abstract_Floor_Penetration
use classes_visitable_walls, only: Abstract_Visitable_Walls
use procedures_walls_factory, only: walls_create => create, walls_destroy => destroy
use classes_min_distance, only: Abstract_Min_Distance
use procedures_hard_core_factory, only: hard_core_create => create, hard_core_destroy => destroy
use types_raw_coordinates, only: Concrete_Raw_Coordinates
use types_component_coordinates_reader_selector, only: Component_Coordinates_Reader_Selector
use procedures_complete_coordinates_reader, only: complete_coordinates_read, &
    complete_coordinates_deallocate
use procedures_property_inquirers, only: use_walls, property_num_components => num_components, &
    property_i_component => i_component
use procedures_plmc_iterations, only: plmc_set_num_snaps
use procedures_plmc_help, only: plmc_catch_density_help

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
    class(Abstract_Visitable_Walls), allocatable :: visitable_walls
    class(Abstract_Floor_Penetration), allocatable :: floor_penetration
    class(Abstract_Min_Distance), allocatable :: wall_min_distance
    real(DP) :: box_size(num_dimensions)
    integer :: num_components, i_component, i_particle
    type(Component_Coordinates_Reader_Selector) :: selector
    integer :: num_snaps, i_snap
    character(len=:), allocatable :: snap_filename
    real(DP), allocatable :: nums_inside(:)
    real(DP) :: avg_num, rms_num
    type(Concrete_Raw_Coordinates) :: raw_coordinates
    type(json_file) :: generating_data, exploring_data

    call plmc_catch_density_help()
    call json_data_create_input(generating_data, 1)
    call plmc_set_num_snaps(num_snaps, generating_data)
    call box_create(periodic_box, generating_data, environment_prefix)
    call walls_create(floor_penetration, generating_data, environment_prefix)
    call hard_core_create(wall_min_distance, use_walls(floor_penetration), generating_data, &
        environment_prefix//"Walls.")
    call walls_create(visitable_walls, periodic_box, floor_penetration, wall_min_distance, &
        generating_data, environment_prefix)
    num_components = property_num_components(generating_data, mixture_prefix)
    call json_data_destroy_input(generating_data)

    call json_data_create_input(exploring_data, 2)
    call box_create(parallelepiped_domain, periodic_box, visitable_walls, .true., exploring_data, &
        "Density.")
    i_component = property_i_component(exploring_data, density_prefix)
    call check_in_range("density", num_components, "i_component", i_component)
    call json_data_destroy_input(exploring_data)

    selector%read_positions = .true.
    selector%read_orientations = .false.
    allocate(nums_inside(num_snaps))
    nums_inside = 0._DP
    do i_snap = 1, num_snaps
        call create_filename_from_argument(snap_filename, i_snap + 2)
        call complete_coordinates_read(box_size, raw_coordinates, num_components, i_component, &
            snap_filename, selector)
        call periodic_box%set(box_size)
        do i_particle = 1, size(raw_coordinates%positions, 2)
            if (parallelepiped_domain%is_inside(raw_coordinates%positions(:, i_particle))) then
                nums_inside(i_snap) = nums_inside(i_snap) + 1._DP
            end if
        end do
        call complete_coordinates_deallocate(raw_coordinates)
        if (allocated(snap_filename)) deallocate(snap_filename)
    end do

    avg_num = sum(nums_inside) / size(nums_inside)
    write(output_unit, *) "domain.density.avg", avg_num / product(parallelepiped_domain%get_size())
    rms_num = sqrt(sum((nums_inside - avg_num)**2) / (size(nums_inside) - 1))
    write(output_unit, *) "domain.density.rms", rms_num / product(parallelepiped_domain%get_size())
    deallocate(nums_inside)

    call box_destroy(parallelepiped_domain)
    call walls_destroy(visitable_walls)
    call hard_core_destroy(wall_min_distance)
    call walls_destroy(floor_penetration)
    call box_destroy(periodic_box)

end program density
