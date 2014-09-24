!> \brief Calculate and print the distribution function between 2 types.

!> \f[
!>      g_{ij}(r) = \frac{\braket{N_\text{pair}(r_{ij}) / N_i}}
!>                       {\braket{\rho_j} \DeltaV_\text{sphere}(r_{ij})}
!> \f]

program between_radial_distribution

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit, error_unit
use data_box, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: data_filename, data_post_filename, report_filename, &
                       test_file_exists, test_data_found
use module_geometry, only: set_geometry
use module_physics_micro, only: sphere_volume, PBC_distance
use module_arguments, only: arg_to_file

implicit none
   
    logical :: take_snapshot
    real(DP), dimension(:), allocatable :: Box_size

    character(len=4096) :: type1_name, type2_name
    integer :: type1_num_particles, type2_num_particles
    integer :: type1_snap_factor, type2_snap_factor
    real(DP) :: type1_density, type2_density
    integer :: type1_positions_unit, type2_positions_unit
    
    integer :: num_distribution, i_distribution
    integer :: num_thermalisation_steps
    integer :: num_equilibrium_steps, i_step, num_common_steps
    integer :: type1_i_particle, type2_i_particle
    real(DP) :: distance_12, distance_max, delta
    real(DP) :: distance_i, distance_minus, distance_plus
    real(DP), dimension(:), allocatable :: distribution_step
    real(DP), dimension(:), allocatable :: distribution_function
    real(DP), dimension(:, :), allocatable :: type1_positions, type2_positions
    
    type(json_file) :: data_json, data_post_json, report_json
    character(len=4096) :: data_name
    logical :: found
    character(len=4096) :: filename
    integer :: length
    integer :: report_unit, distrib_unit
    real(DP) :: time_start, time_end

    character(len=:), allocatable :: geometry
    
    call json_initialize()
    
    call test_file_exists(data_filename)
    call data_json%load_file(filename = data_filename)
    
    data_name = "Distribution.take snapshot"
    call data_json%get(data_name, take_snapshot, found)
    call test_data_found(data_name, found)
    
    if (.not.take_snapshot) stop "No snap shots taken."
    
    call test_file_exists(report_filename)
    call report_json%load_file(filename = report_filename)

    data_name = "System.Box.geometry"
    call report_json%get(data_name, geometry, found)
    call test_data_found(data_name, found)
    call set_geometry(geometry)
    if (allocated(geometry)) deallocate(geometry)

    call report_json%destroy()
    
    data_name = "Box.size"
    call data_json%get(data_name, Box_size, found)
    call test_data_found(data_name, found)
    if (size(Box_size) /= num_dimensions) error stop "Box size dimension"
    
    data_name = "Monte Carlo.number of thermalisation steps"
    call data_json%get(data_name, num_thermalisation_steps, found)
    call test_data_found(data_name, found)
    
    data_name = "Monte Carlo.number of equilibrium steps"
    call data_json%get(data_name, num_equilibrium_steps, found)
    call test_data_found(data_name, found)

    call data_json%destroy()

    call test_file_exists(data_post_filename)
    call data_post_json%load_file(filename = data_post_filename)
    
    data_name = "Distribution.delta"
    call data_post_json%get(data_name, delta, found)
    call test_data_found(data_name, found)

    call data_post_json%destroy()
    
    distance_max = norm2(Box_size/2._DP)
    num_distribution = int(distance_max/delta)
    allocate(distribution_step(num_distribution))
    allocate(distribution_function(num_distribution))
    
    call arg_to_file(1, filename, length)
    open(newunit=type1_positions_unit, recl=4096, file=filename(1:length), status='old', action='read')
    read(type1_positions_unit, *) type1_name, type1_num_particles, type1_snap_factor
    write(output_unit, *) "type 1: ", trim(type1_name), type1_num_particles, type1_snap_factor
    allocate(type1_positions(num_dimensions, type1_num_particles))
    type1_density = type1_num_particles / product(Box_size)
    
    call arg_to_file(2, filename, length)
    open(newunit=type2_positions_unit, recl=4096, file=filename(1:length), status='old', action='read')
    read(type2_positions_unit, *) type2_name, type2_num_particles, type2_snap_factor
    write(output_unit, *) "type 2: ", trim(type2_name), type2_num_particles, type2_snap_factor
    allocate(type2_positions(num_dimensions, type2_num_particles))
    type2_density = type2_num_particles / product(Box_size)
    
    write(output_unit, *) "Start !"
    call cpu_time(time_start)
    distribution_function(:) = 0._DP
    num_common_steps = 0
    do i_step = num_thermalisation_steps + 1, num_thermalisation_steps + num_equilibrium_steps
        
        if (modulo(i_step, type1_snap_factor) == 0) then
            do type1_i_particle = 1, type1_num_particles
                read(type1_positions_unit, *) type1_positions(:, type1_i_particle)
            end do
        end if
        
        if (modulo(i_step, type2_snap_factor) == 0) then
            do type2_i_particle = 1, type2_num_particles
                read(type2_positions_unit, *) type2_positions(:, type2_i_particle)
            end do
        end if
        
        if (modulo(i_step, type1_snap_factor) == 0 .and. modulo(i_step, type2_snap_factor) == 0) then
            num_common_steps = num_common_steps + 1
            distribution_step(:) = 0._DP
            do type1_i_particle = 1, type1_num_particles
                do type2_i_particle = 1, type2_num_particles
                    distance_12 = PBC_distance(Box_size, &
                                               type1_positions(:, type1_i_particle), &
                                               type2_positions(:, type2_i_particle))
                    i_distribution = int(distance_12/delta)
                    distribution_step(i_distribution) = distribution_step(i_distribution) + 1._DP
                end do
            end do
            distribution_function(:) = distribution_function(:) + distribution_step(:)
        end if
    
    end do
    call cpu_time(time_end)
    write(output_unit, *) "Finish !"
    
    deallocate(type2_positions)
    close(type2_positions_unit)
    
    deallocate(type1_positions)
    close(type1_positions_unit)
    
    open(newunit=distrib_unit, file=trim(type1_name)//"-"//trim(type2_name)//"_radial_distribution_function.out", &
         action="write")
    
        distribution_function(:) = distribution_function(:) / real(num_common_steps, DP) / &
                                   real(type1_num_particles, DP)
    
        do i_distribution = 1, num_distribution
        
            distance_i = (real(i_distribution, DP) + 0.5_DP) * delta
            distance_minus = real(i_distribution, DP) * delta
            distance_plus = real(i_distribution + 1, DP) * delta
            
            distribution_function(i_distribution) = distribution_function(i_distribution) / &
                type2_density / (sphere_volume(distance_plus) - sphere_volume(distance_minus))
            write(distrib_unit, *) distance_i, distribution_function(i_distribution)
            
        end do
        
    close(distrib_unit)

    open(newunit=report_unit, file=trim(type1_name)//"-"//trim(type2_name)//"_radial_distribution_report.txt", &
         action="write")
         write(report_unit, *) "Duration =", (time_end - time_start) / 60._DP, "min"
    close(report_unit)
    
    deallocate(distribution_function)
    deallocate(distribution_step)

end program between_radial_distribution
