program plmc_domain_z_distribution

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_constants, only: num_dimensions
use data_strings, only: max_line_length
use data_prefixes, only: environment_prefix
use json_module, only: json_file
use procedures_errors, only: error_exit, warning_continue
use procedures_checks, only: check_data_found, check_string_not_empty
use procedures_command_arguments, only: create_filename_from_argument
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_box_factory, only: box_create => create, box_destroy => destroy
use procedures_coordinates_reader, only: create_positions_from_file
use procedures_plmc_factory, only: plmc_create, plmc_destroy
use procedures_property_inquirers, only: periodicity_is_xy

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
    integer :: i_particle
    logical :: coordinates_written
    integer :: num_snaps, i_snap
    character(len=:), allocatable :: snap_filename
    real(DP), dimension(num_dimensions) :: domain_origin, domain_size
    real(DP), dimension(:, :), allocatable :: positions
    real(DP) :: delta
    integer :: bins_bounds(2)
    real(DP), allocatable :: bins_function(:)
    integer :: i_z
    integer :: bins_unit
    character(len=:), allocatable :: bins_filename

    type(json_file) :: input_data, post_data
    character(len=:), allocatable :: data_field
    logical :: data_found

    call plmc_create(input_data, command_argument_count() - 1)
    data_field = "Output.Coordinates.write"
    call input_data%get(data_field, coordinates_written, data_found)
    call check_data_found(data_field, data_found)
    if (.not.coordinates_written) call error_exit("Coordinates weren't written.")
    num_snaps = command_argument_count() - 2
    call box_create(periodic_box, input_data, environment_prefix)
    call plmc_destroy(input_data)

    call plmc_create(post_data, command_argument_count())
    call box_create(parallelepiped_domain, periodic_box, .true., post_data, "Z Distribution.")
    if (.not.periodicity_is_xy(periodic_box)) then
        call warning_continue("Periodicity is not XY.")
    end if
    data_field = "Z Distribution.delta"
    call post_data%get(data_field, delta, data_found)
    call check_data_found(data_field, data_found)
    data_field = "Z Distribution.file name"
    call post_data%get(data_field, bins_filename, data_found)
    call check_data_found(data_field, data_found)
    call check_string_not_empty(data_field, bins_filename)
    call plmc_destroy(post_data)
    deallocate(data_field)

    domain_origin = parallelepiped_domain%get_origin()
    domain_size = parallelepiped_domain%get_size()
    bins_bounds = [nint((domain_origin(3) - domain_size(3)/2)/delta), &
                   nint((domain_origin(3) + domain_size(3)/2)/delta)]
    allocate(bins_function(bins_bounds(1):bins_bounds(2)))

    bins_function = 0._DP
    do i_snap = 1, num_snaps
        call create_filename_from_argument(snap_filename, i_snap)
        call create_positions_from_file(positions, snap_filename)
        do i_particle = 1, size(positions, 2)
            if (parallelepiped_domain%is_inside(positions(:, i_particle))) then
                i_z = nint(positions(3, i_particle)/delta)
                bins_function(i_z) = bins_function(i_z) + 1._DP
            end if
        end do
        deallocate(positions)
        deallocate(snap_filename)
    end do

    bins_function = bins_function / real(num_snaps, DP) / product(domain_size(1:2)) / delta

    open(unit=bins_unit, recl=max_line_length, file=bins_filename, action="write")
    deallocate(bins_filename)
    write(bins_unit, *) "#  z   distribution"
    do i_z = lbound(bins_function, 1), ubound(bins_function, 1)
        write(bins_unit, *) real(i_z, DP) * delta, bins_function(i_z)
    end do
    close(bins_unit)

    deallocate(bins_function)
    call box_destroy(periodic_box)
    call box_destroy(parallelepiped_domain)

end program plmc_domain_z_distribution
