!Calculate the radial distribution function from snap shots of 1 component
!> \[
!>      g(r) = \frac{\mathrm{d}N}{\mathrm{d}r} \frac{1}{\rho S(r)}
!> \]
!> with \( S(r) \) the area of a sphere.

program radial_intra

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_strings, only: max_line_length
use data_input_prefixes, only: environment_prefix
use json_module, only: json_file
use procedures_json_data_factory, only: json_data_create_input => create_input, &
    json_data_destroy_input => destroy_input
use procedures_errors, only: error_exit, warning_continue
use procedures_checks, only: check_data_found, check_positive, check_string_not_empty
use procedures_geometry, only: sphere_surface
use procedures_command_arguments, only: create_filename_from_argument
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_box_factory, only: box_create => create, box_destroy => destroy
use procedures_coordinates_reader, only: create_positions_from_file
use procedures_property_inquirers, only: periodicity_is_xyz
use types_radial_distribution_component, only: Concrete_Radial_Distribution_Component
use procedures_plmc_help, only: plmc_catch_radial_help

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    type(Concrete_Radial_Distribution_Component) :: component
    integer :: i_particle, j_particle
    logical :: coordinates_written
    integer :: num_snaps, i_snap
    character(len=:), allocatable :: snap_filename
    real(DP) :: max_distance, delta_distance, distance_ij, distance_i
    real(DP), allocatable :: bins_snap(:), bins_function(:)
    integer :: i_bin
    integer :: bins_unit
    character(len=:), allocatable :: bins_filename

    type(json_file) :: generating_data, exploring_data
    character(len=:), allocatable :: data_field
    logical :: data_found

    call plmc_catch_radial_help()
    num_snaps = command_argument_count() - 2
    if (num_snaps == 0) call error_exit("No snaps given.")

    call json_data_create_input(generating_data, 1)
    data_field = "Output.Coordinates.write"
    call generating_data%get(data_field, coordinates_written, data_found)
    call check_data_found(data_field, data_found)
    if (.not.coordinates_written) call error_exit("Coordinates weren't written.")
    call box_create(periodic_box, generating_data, environment_prefix)
    if (.not.periodicity_is_xyz(periodic_box)) then
        call warning_continue("Periodicity is not XYZ.")
    end if
    call json_data_destroy_input(generating_data)
    max_distance = minval(periodic_box%get_size())/2._DP

    call json_data_create_input(exploring_data, 2)
    data_field = "Radial.delta"
    call exploring_data%get(data_field, delta_distance, data_found)
    call check_data_found(data_field, data_found)
    call check_positive("radial_distribution", "delta_distance", delta_distance)
    data_field = "Radial.file name"
    call exploring_data%get(data_field, bins_filename, data_found)
    call check_data_found(data_field, data_found)
    call check_string_not_empty(data_field, bins_filename)
    deallocate(data_field)
    call json_data_destroy_input(exploring_data)

    allocate(bins_snap(nint(max_distance/delta_distance)))
    allocate(bins_function(size(bins_snap)))

    component%num_particles_sum = 0
    bins_function = 0._DP
    do i_snap = 1, num_snaps
        call create_filename_from_argument(snap_filename, i_snap + 2)
        call create_positions_from_file(component%positions, snap_filename)
        component%num_particles = size(component%positions, 2)
        bins_snap = 0._DP
        do i_particle = 1, component%num_particles
            do j_particle = i_particle + 1, component%num_particles
                distance_ij = periodic_box%distance(component%positions(:, i_particle), &
                    component%positions(:, j_particle))
                if (distance_ij < max_distance) then
                    i_bin = nint(distance_ij/delta_distance)
                    bins_snap(i_bin) = bins_snap(i_bin) + 1._DP
                end if
            end do
        end do
        deallocate(component%positions)
        deallocate(snap_filename)
        if (component%num_particles > 0) then
            bins_function = bins_function + bins_snap/real(component%num_particles, DP)
            component%num_particles_sum = component%num_particles_sum + component%num_particles
        end if
    end do
    deallocate(bins_snap)

    bins_function = 2._DP * bins_function / real(num_snaps, DP)
    component%density = real(component%num_particles_sum, DP) / real(num_snaps) / &
        product(periodic_box%get_size())
    call box_destroy(periodic_box)

    open(newunit=bins_unit, recl=max_line_length, file=bins_filename, action="write")
    deallocate(bins_filename)
    do i_bin = 1, size(bins_function)
        distance_i = real(i_bin, DP) * delta_distance
        bins_function(i_bin) = bins_function(i_bin) / delta_distance / component%density / &
            sphere_surface(distance_i)
        write(bins_unit, *) distance_i, bins_function(i_bin)
    end do
    close(bins_unit)
    deallocate(bins_function)

end program radial_intra
