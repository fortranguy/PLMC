!Calculate the radial distribution function from snap shots
!> \[
!>      g(r) = \frac{\mathrm{d}N}{\mathrm{d}r} \frac{1}{\rho S(r)}
!> \]
!> with \( S(r) \) the area of a sphere.


program plmc_radial_distribution

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_constants, only: max_line_length
use data_wrappers_prefix, only: environment_prefix
use json_module, only: json_file, json_initialize
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found, check_positive, check_string_not_empty
use procedures_geometry, only: sphere_surface
use procedures_command_arguments, only: create_filename_from_argument
use class_periodic_box, only: Abstract_Periodic_Box
use procedures_environment_factory, only: environment_create, environment_destroy
use procedures_coordinates_micro, only: create_coordinates_from_file
use procedures_plmc_factory, only: plmc_create, plmc_destroy

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    integer :: num_particles_sum, num_particles_snap, i_particle, j_particle
    real(DP) :: density_average
    logical :: coordinates_written
    real(DP), allocatable :: positions(:, :)
    real(DP) :: max_distance, delta_distance, distance_ij, distance_i
    integer :: num_snaps, i_snap
    character(len=:), allocatable :: snap_filename
    integer :: num_bins, i_bin
    real(DP), allocatable :: bins_snap(:), bins_function(:)
    integer :: bins_unit
    character(len=:), allocatable :: bins_filename

    type(json_file) :: input_data, post_data
    character(len=:), allocatable :: data_field
    logical :: data_found

    call json_initialize()

    call plmc_create(input_data, command_argument_count() - 1)
    data_field = "Output.Coordinates.write"
    call input_data%get(data_field, coordinates_written, data_found)
    call check_data_found(data_field, data_found)
    if (.not.coordinates_written) then
        call error_exit("Coordinates weren't written.")
    end if
    num_snaps = command_argument_count() - 2
    call environment_create(periodic_box, input_data, environment_prefix)!warning: if not 3D
    call plmc_destroy(input_data)
    max_distance = minval(periodic_box%get_size())/2._DP

    call plmc_create(post_data, command_argument_count())
    data_field = "Distribution.delta"
    call post_data%get(data_field, delta_distance, data_found)
    call check_data_found(data_field, data_found)
    call check_positive("radial_distribution", "delta_distance", delta_distance)
    data_field = "Distribution.file name"
    call post_data%get(data_field, bins_filename, data_found)
    call check_data_found(data_field, data_found)
    call check_string_not_empty(data_field, bins_filename)
    call plmc_destroy(post_data)
    deallocate(data_field)

    num_bins = int(max_distance/delta_distance) + 1
    allocate(bins_snap(num_bins))
    allocate(bins_function(num_bins))

    num_particles_sum = 0
    bins_function = 0._DP
    do i_snap = 1, num_snaps
        call create_filename_from_argument(snap_filename, i_snap)
        call create_coordinates_from_file(positions, snap_filename)
        num_particles_snap = size(positions, 2)
        bins_snap = 0._DP
        do i_particle = 1, num_particles_snap
            do j_particle = i_particle + 1, num_particles_snap
                distance_ij = periodic_box%distance(positions(:, i_particle), &
                    positions(:, j_particle))
                if (distance_ij < max_distance) then
                    i_bin = int(distance_ij/delta_distance) + 1
                    bins_snap(i_bin) = bins_snap(i_bin) + 1._DP
                end if
            end do
        end do
        deallocate(positions)
        deallocate(snap_filename)
        if (num_particles_snap > 0) then
            bins_function = bins_function + bins_snap/real(num_particles_snap, DP)
            num_particles_sum = num_particles_sum + num_particles_snap
        end if
    end do
    deallocate(bins_snap)

    bins_function = 2._DP * bins_function / real(num_snaps, DP)
    density_average = real(num_particles_sum, DP) / real(num_snaps) / &
        product(periodic_box%get_size())
    call environment_destroy(periodic_box)

    open(newunit=bins_unit, recl=max_line_length, file=bins_filename, action="write")
    do i_bin = 1, num_bins
        distance_i = (real(i_bin, DP) - 0.5) * delta_distance
        bins_function(i_bin) = bins_function(i_bin) / delta_distance / density_average / &
            sphere_surface(distance_i)
        write(bins_unit, *) distance_i, bins_function(i_bin)
    end do
    close(bins_unit)
    deallocate(bins_function)

end program plmc_radial_distribution
