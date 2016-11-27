program radial

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use data_input_prefixes, only: environment_prefix, mixture_prefix
use data_arguments, only: i_generating, i_exploring, num_json_arguments
use json_module, only: json_file
use procedures_json_data_factory, only: json_data_create => create, json_data_destroy => destroy
use procedures_errors, only: warning_continue
use procedures_command_arguments, only: create_filename_from_argument
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_boxes_factory, only: boxes_create => create, boxes_destroy => destroy
use procedures_environment_inquirers, only: periodicity_is_xyz
use procedures_mixture_inquirers, only: property_num_components => num_components
use procedures_complete_coordinates_reader, only: complete_coordinates_read => read
use procedures_plmc_iterations, only: plmc_set_num_snaps
use classes_radial_explorer, only: Abstract_Radial_Explorer
use procedures_radial_explorer_factory, only: radial_explorer_create => create, &
    radial_explorer_destroy => destroy
use procedures_plmc_help, only: plmc_catch_radial_help

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_boxes(:)
    real(DP), dimension(num_dimensions) :: max_box_size, box_size
    class(Abstract_Radial_Explorer), allocatable :: radial_explorer
    integer :: num_components, num_snaps, i_snap
    character(len=:), allocatable :: snap_filename
    type(json_file) :: generating_data, exploring_data

    call plmc_catch_radial_help()
    call json_data_create(generating_data, i_generating)
    call plmc_set_num_snaps(num_snaps, generating_data)
    call boxes_create(periodic_boxes, generating_data, environment_prefix, unique=.true.)
    if (.not.periodicity_is_xyz(periodic_boxes(1))) then
        call warning_continue("Periodicity is not XYZ.")
    end if
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
    call radial_explorer_create(radial_explorer, max_box_size, num_components, exploring_data)
    call json_data_destroy(exploring_data)

    do i_snap = 1, num_snaps
        call create_filename_from_argument(snap_filename, num_json_arguments + i_snap)
        call radial_explorer%read_and_fill(periodic_boxes(1), snap_filename)
        if (allocated(snap_filename)) deallocate(snap_filename)
    end do
    call radial_explorer%write(num_snaps)

    call radial_explorer_destroy(radial_explorer)
    call boxes_destroy(periodic_boxes)

end program radial
