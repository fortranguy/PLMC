program plmc_domain_density

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_prefixes, only: environment_prefix
use json_module, only: json_file, json_initialize
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use procedures_command_arguments, only: create_filename_from_argument
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_environment_factory, only: environment_create, environment_destroy
use procedures_coordinates_reader, only: create_positions_from_file
use procedures_plmc_factory, only: plmc_create, plmc_destroy

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
    integer :: i_particle
    logical :: coordinates_written
    integer :: num_snaps, i_snap
    character(len=:), allocatable :: snap_filename
    integer :: num_inside_sum
    real(DP), dimension(:, :), allocatable :: positions

    type(json_file) :: input_data, post_data
    character(len=:), allocatable :: data_field
    logical :: data_found

    call json_initialize()

    call plmc_create(input_data, command_argument_count() - 1)
    data_field = "Output.Coordinates.write"
    call input_data%get(data_field, coordinates_written, data_found)
    call check_data_found(data_field, data_found)
    deallocate(data_field)
    if (.not.coordinates_written) call error_exit("Coordinates weren't written.")
    num_snaps = command_argument_count() - 2
    call environment_create(periodic_box, input_data, environment_prefix)
    call plmc_destroy(input_data)

    call plmc_create(post_data, command_argument_count())
    call environment_create(parallelepiped_domain, .true., periodic_box, post_data, "Density.")
    call plmc_destroy(post_data)

    num_inside_sum = 0
    do i_snap = 1, num_snaps
        call create_filename_from_argument(snap_filename, i_snap)
        call create_positions_from_file(positions, snap_filename)
        do i_particle = 1, size(positions, 2)
            if (parallelepiped_domain%is_inside(positions(:, i_particle))) then
                num_inside_sum = num_inside_sum + 1
            end if
        end do
        deallocate(positions)
        deallocate(snap_filename)
    end do

    write(output_unit, *) "domain.density", real(num_inside_sum, DP) / real(num_snaps, DP) / &
        product(parallelepiped_domain%get_size())
    call environment_destroy(periodic_box)
    call environment_destroy(parallelepiped_domain)

end program plmc_domain_density
