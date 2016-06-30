program density

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_input_prefixes, only: environment_prefix
use json_module, only: json_file
use procedures_json_data_factory, only: json_data_create_input => create_input, &
    json_data_destroy_input => destroy_input
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use procedures_command_arguments, only: create_filename_from_argument
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_box_factory, only: box_create => create, box_destroy => destroy
use classes_floor_penetration, only: Abstract_Floor_Penetration
use classes_visitable_walls, only: Abstract_Visitable_Walls
use procedures_walls_factory, only: walls_create => create, walls_destroy => destroy
use procedures_coordinates_reader, only: create_positions_from_file
use procedures_plmc_help, only: plmc_catch_density_help

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
    class(Abstract_Visitable_Walls), allocatable :: visitable_walls
    class(Abstract_Floor_Penetration), allocatable :: floor_penetration
    integer :: i_particle
    logical :: coordinates_written
    integer :: num_snaps, i_snap
    character(len=:), allocatable :: snap_filename
    real(DP), allocatable :: nums_inside(:)
    real(DP) :: avg_num, rms_num
    real(DP), dimension(:, :), allocatable :: positions

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
    deallocate(data_field)
    if (.not.coordinates_written) call error_exit("Coordinates weren't written.")
    call box_create(periodic_box, generating_data, environment_prefix)
    call walls_create(floor_penetration, generating_data, environment_prefix)
    call walls_create(visitable_walls, periodic_box, floor_penetration, generating_data, &
        environment_prefix)
    call json_data_destroy_input(generating_data)
    call json_data_create_input(exploring_data, 2)
    call box_create(parallelepiped_domain, periodic_box, visitable_walls, .true., exploring_data, &
        "Density.")
    call json_data_destroy_input(exploring_data)

    allocate(nums_inside(num_snaps))
    nums_inside = 0._DP
    do i_snap = 1, num_snaps
        call create_filename_from_argument(snap_filename, i_snap + 2)
        call create_positions_from_file(positions, snap_filename)
        do i_particle = 1, size(positions, 2)
            if (parallelepiped_domain%is_inside(positions(:, i_particle))) then
                nums_inside(i_snap) = nums_inside(i_snap) + 1._DP
            end if
        end do
        deallocate(positions)
        deallocate(snap_filename)
    end do

    avg_num = sum(nums_inside) / size(nums_inside)
    write(output_unit, *) "domain.density.avg", avg_num / product(parallelepiped_domain%get_size())
    rms_num = sqrt(sum((nums_inside - avg_num)**2) / (size(nums_inside) - 1))
    write(output_unit, *) "domain.density.rms", rms_num / product(parallelepiped_domain%get_size())
    deallocate(nums_inside)

    call box_destroy(parallelepiped_domain)
    call walls_destroy(visitable_walls)
    call walls_destroy(floor_penetration)
    call box_destroy(periodic_box)

end program density
