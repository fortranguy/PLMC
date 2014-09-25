!> \brief Calculate and write the distribution function

!> \f[
!>      g_\parallel(r^{2D}) = \frac{\braket{N_\text{pair}} S}
!>                                 {\braket{N}_{\Delta z} 2\pi r^{2D} \Delta r^{2D}}
!> \f]

program parallel_distribution

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit, error_unit
use data_precisions, only: real_zero
use data_constants, only: PI
use data_box, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: data_filename, data_post_filename, report_filename, &
                       test_file_exists, test_data_found
use module_geometry, only: set_geometry, geometry
use module_physics_micro, only: sphere_volume, PBC_vector
use module_arguments, only: arg_to_file

implicit none

    logical :: take_snapshot
    character(len=:), allocatable :: Box_geometry
    real(DP), dimension(:), allocatable :: Box_size
    real(DP) :: Box_height

    character(len=4096) :: name
    integer :: num_particles, num_particles_inside_sum, num_particles_inside
    real(DP) :: num_particles_inside_avg
    integer :: snap_factor
    real(DP) :: diameter
    integer :: positions_unit, distrib_unit
    
    integer :: num_thermalisation_steps
    integer :: num_equilibrium_steps, i_step, num_steps
    integer :: i_particle, j_particle
    real(DP) :: distance_max, delta
    real(DP), dimension(num_dimensions) :: vector_ij
    real(DP) :: distance_i_distribution
    integer :: num_distribution, i_distribution
    real(DP), dimension(:), allocatable :: distribution_step
    real(DP), dimension(:), allocatable :: distribution_function
    real(DP), dimension(:, :), allocatable :: positions
    real(DP) :: height_delta
    real(DP) :: height_ratio, height, z_min, z_max
    
    type(json_file) :: data_json, data_post_json, report_json
    character(len=4096) :: data_name
    logical :: found
    character(len=4096) :: filename
    integer :: length
    integer :: report_unit
    real(DP) :: time_start, time_end
    
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
    call report_json%get(data_name, Box_geometry, found)
    call test_data_found(data_name, found)
    call set_geometry(Box_geometry)
    if (allocated(Box_geometry)) deallocate(Box_geometry)
    
    if (.not. geometry%slab) error stop "For slab geometry only."

    call report_json%destroy()
    
    data_name = "Box.size"
    call data_json%get(data_name, Box_size, found)
    call test_data_found(data_name, found)
    if (size(Box_size) /= 2) error stop "Box size dimension"
    
    data_name = "Box.height"
    call data_json%get(data_name, Box_height, found)
    call test_data_found(data_name, found)
    
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
    
    data_name = "Distribution.Parallel.height ratio"
    call data_post_json%get(data_name, height_ratio, found)
    call test_data_found(data_name, found)
    
    data_name = "Distribution.Parallel.height delta"
    call data_post_json%get(data_name, height_delta, found)
    call test_data_found(data_name, found)

    call data_post_json%destroy()
    
    distance_max = norm2(Box_size(1:2)/2._DP)
    num_distribution = int(distance_max/delta) + 1
    allocate(distribution_step(num_distribution))
    allocate(distribution_function(num_distribution))
    
    call arg_to_file(1, filename, length)
    open(newunit=positions_unit, recl=4096, file=filename(1:length), status='old', action='read')
    
    read(positions_unit, *) name, num_particles, snap_factor
    write(output_unit, *) trim(name), num_particles, snap_factor
    diameter = 1._DP
    
    allocate(positions(num_dimensions, num_particles))
    
    height = height_ratio * Box_height ! 0. <= height_ratio <= 1.
    
    if (diameter/2._DP <= height .and. height < height_delta + diameter/2._DP) then
        z_min = height
        z_max = z_min + height_delta
    else if (abs(height_ratio - 0._DP) < real_zero) then
        z_min = diameter/2._DP
        z_max = z_min + height_delta
    else if (Box_height - (height_delta + diameter/2._DP) < height .and. &
             height <= Box_height - diameter/2._DP) then
        z_max = height
        z_min = z_max - height_delta
    else if (abs(height_ratio - 1._DP) < real_zero) then
        z_max = Box_height - diameter/2._DP
        z_min = z_max - height_delta
    else
        z_min = height - height_delta/2._DP
        z_max = height + height_delta/2._DP
    end if

    write(output_unit, *) "Start !"
    call cpu_time(time_start)
    num_particles_inside_sum = 0
    distribution_function(:) = 0._DP
    num_steps = 0
    do i_step = num_thermalisation_steps + 1, num_thermalisation_steps + num_equilibrium_steps
        if (modulo(i_step, snap_factor) == 0) then
        
            num_steps = num_steps + 1
        
            do i_particle = 1, num_particles
                read(positions_unit, *) positions(:, i_particle)
            end do

            ! Fill
            num_particles_inside = 0
            distribution_step(:) = 0
            do i_particle = 1, num_particles
                if (z_min < positions(3, i_particle) .and. positions(3, i_particle) < z_max) then
                    num_particles_inside = num_particles_inside + 1
                    do j_particle = i_particle + 1, num_particles
                        if (z_min < positions(3, j_particle) .and. positions(3, j_particle) < z_max) then
                            vector_ij = PBC_vector(Box_size, positions(:, i_particle), &
                                                             positions(:, j_particle))
                            i_distribution =  int(norm2(vector_ij(1:2)) / delta) + 1
                            distribution_step(i_distribution) = distribution_step(i_distribution) + 1._DP
                        end if
                    end do
                end if
            end do
            
            num_particles_inside_sum = num_particles_inside_sum + num_particles_inside
            distribution_function(:) = distribution_function(:) + distribution_step(:)
        
        end if
    end do
    call cpu_time(time_end)
    write(output_unit, *) "Finish !"

    close(positions_unit)
    deallocate(positions)

    open(newunit=distrib_unit, file=trim(name)//"_parallel_distribution_function.out", &
         action="write")
         
        num_particles_inside_avg = real(num_particles_inside_sum, DP) / real(num_steps, DP)
        distribution_function(:) = 2._DP * distribution_function(:) / real(num_steps, DP)
        distribution_function(:) = distribution_function(:) * product(Box_size(1:2)) / &
                                   num_particles_inside_avg**2 / (2._DP*PI) / delta
    
        do i_distribution = 1, num_distribution
            distance_i_distribution = (real(i_distribution, DP) - 0.5_DP) * delta
            write(distrib_unit, *) distance_i_distribution, distribution_function(i_distribution) / &
                                                            distance_i_distribution
        end do
        
    close(distrib_unit)
    
    open(newunit=report_unit, file=trim(name)//"_parallel_distribution_report.txt", &
         action="write")
        write(report_unit, *) "Duration =", (time_end - time_start) / 60._DP, "min"
    close(report_unit)
    
    deallocate(distribution_function)
    deallocate(distribution_step)
    
    deallocate(Box_size)

end program parallel_distribution
