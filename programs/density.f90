program density

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_input_prefixes, only: environment_prefix, mixture_prefix
use data_constants, only: num_dimensions
use json_module, only: json_file
use procedures_json_data_factory, only: json_data_create_input => create_input, &
    json_data_destroy_input => destroy_input
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found, check_in_range
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
use procedures_complete_coordinates_reader, only: complete_coordinates_read
use procedures_property_inquirers, only: use_walls
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
    logical :: coordinates_written
    integer :: num_snaps, i_snap
    character(len=:), allocatable :: snap_filename
    real(DP), allocatable :: nums_inside(:)
    real(DP) :: avg_num, rms_num
    type(Concrete_Raw_Coordinates) :: raw_coordinates

    type(json_file) :: generating_data, exploring_data
    character(len=:), allocatable :: data_field
    logical :: data_found

    call plmc_catch_density_help()
    num_snaps = command_argument_count() - 2
    if (num_snaps == 0) call error_exit("No snaps given.")

    call json_data_create_input(generating_data, 1)
    data_field = "Output.Coordinates.write"
    call generating_data%get(data_field, coordinates_written, data_found)
    call check_data_found(data_field, data_found)
    if (.not.coordinates_written) call error_exit("Coordinates weren't written.")
    call box_create(periodic_box, generating_data, environment_prefix)
    call walls_create(floor_penetration, generating_data, environment_prefix)
    call hard_core_create(wall_min_distance, use_walls(floor_penetration), generating_data, &
        environment_prefix//"Walls.")
    call walls_create(visitable_walls, periodic_box, floor_penetration, wall_min_distance, &
        generating_data, environment_prefix)
    data_field = mixture_prefix//"number of components"
    call generating_data%get(data_field, num_components, data_found)
    call check_data_found(data_field, data_found)

    call json_data_create_input(exploring_data, 2)
    call box_create(parallelepiped_domain, periodic_box, visitable_walls, .true., exploring_data, &
        "Density.")
    data_field = "Density.component number"
    call exploring_data%get(data_field, i_component, data_found)
    call check_data_found(data_field, data_found)
    call check_in_range("density", num_components, "i_component", i_component)
    call json_data_destroy_input(exploring_data)

    selector%read_positions = .true.
    selector%read_orientations = .false.
    deallocate(data_field)
    call json_data_destroy_input(generating_data)

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
        deallocate(raw_coordinates%orientations)
        deallocate(raw_coordinates%positions)
        deallocate(snap_filename)
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
