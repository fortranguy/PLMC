!> \brief Calculate and write the distribution function

!> \f[
!>      g_{ii}(r) = 2 \frac{\braket{\mathrm{d} N_\text{pair}(r) / N_i}}{\braket{\rho_i}
!>                        \mathrm{d} V_\text{sphere}(r)}
!> \f]

program radial_distribution

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit, error_unit
use data_precisions, only: real_zero
use data_box, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: data_filename, data_post_filename, report_filename, &
                       test_file_exists, test_data_found
use module_geometry, only: set_geometry, geometry
use module_physics_micro, only: sphere_volume, PBC_distance, set_bounds
use module_physics_macro, only: test_particles_inside
use module_arguments, only: arg_to_file

implicit none

    logical :: take_snapshot
    real(DP), dimension(:), allocatable :: Box_size
    real(DP) :: Box_height
    real(DP), dimension(:), allocatable :: domain_ratio
    real(DP), dimension(num_dimensions) :: Box_lower_bounds, Box_upper_bounds
    real(DP) :: Volume_inside

    character(len=4096) :: name
    integer :: num_particles
    integer :: snap_factor
    real(DP) :: num_particles_inside
    integer :: num_particles_sum, num_particles_step
    real(DP) :: density_inside
    logical, dimension(:), allocatable :: particles_inside
    integer :: positions_unit
    
    integer :: num_distribution, i_distribution
    integer :: num_thermalisation_steps
    integer :: num_equilibrium_steps, i_step, num_steps
    integer :: i_particle, j_particle
    real(DP) :: distance_ij, distance_max, delta
    real(DP) :: distance_i, distance_i_minus, distance_i_plus
    real(DP) :: cutoff
    real(DP), dimension(:), allocatable :: distribution_step
    real(DP), dimension(:), allocatable :: distribution_function
    real(DP), dimension(:, :), allocatable :: positions
    
    type(json_file) :: data_json, data_post_json, report_json
    character(len=4096) :: data_name
    logical :: found
    character(len=4096) :: filename
    integer :: length
    integer :: report_unit, distrib_unit
    real(DP) :: time_start, time_end

    character(len=:), allocatable :: Box_geometry
    
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

    call report_json%destroy()
    
    data_name = "Box.size"
    call data_json%get(data_name, Box_size, found)
    call test_data_found(data_name, found)
    
    if (geometry%bulk) then
        if (size(Box_size) /= num_dimensions) error stop "Box size dimension"
    else if (geometry%slab) then
        if (size(Box_size) /= 2) error stop "Box size dimension"
        data_name = "Box.height"
        call data_json%get(data_name, Box_height, found)
        call test_data_found(data_name, found)
    end if
    
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

    data_name = "Distribution.Radial.domain ratio"
    call data_post_json%get(data_name, domain_ratio, found)
    call test_data_found(data_name, found)
    if (size(domain_ratio) /= num_dimensions) error stop "domain ratio dimension"
    
    call set_bounds(Box_size, Box_height, domain_ratio, Box_lower_bounds, Box_upper_bounds)
    Volume_inside = product(Box_upper_bounds - Box_lower_bounds)
    
    data_name = "Distribution.Radial.cut off"
    call data_post_json%get(data_name, cutoff, found)
    call test_data_found(data_name, found)

    call data_post_json%destroy()
    
    if (geometry%bulk) then
        distance_max = norm2(Box_size/2._DP)
    else if (geometry%slab) then
        distance_max = min(norm2(Box_size(1:2)), Box_height) /2._DP ! incorrect?
    end if
    
    if (cutoff < real_zero .or. distance_max < cutoff) then
        cutoff = distance_max
    else 
        distance_max = cutoff
    end if
    
    num_distribution = int(distance_max/delta) + 1
    allocate(distribution_step(num_distribution))
    allocate(distribution_function(num_distribution))
    
    call arg_to_file(1, filename, length)
    open(newunit=positions_unit, recl=4096, file=filename(1:length), status='old', action='read')
    
    read(positions_unit, *) name, num_particles, snap_factor
    write(output_unit, *) trim(name), num_particles, snap_factor
    
    allocate(positions(num_dimensions, num_particles))
    allocate(particles_inside(num_particles))

    write(output_unit, *) "Start !"
    call cpu_time(time_start)
    num_particles_sum = 0
    distribution_function(:) = 0._DP
    num_steps = 0    
    do i_step = num_thermalisation_steps + 1, num_thermalisation_steps + num_equilibrium_steps
        
        if (modulo(i_step, snap_factor) == 0) then
        
            num_steps = num_steps + 1
        
            do i_particle = 1, num_particles
                read(positions_unit, *) positions(:, i_particle)
            end do

            call test_particles_inside(Box_lower_bounds, Box_upper_bounds, &
                                       num_particles, positions, particles_inside, &
                                       num_particles_step)
            num_particles_sum = num_particles_sum + num_particles_step

            distribution_step(:) = 0
            do i_particle = 1, num_particles
                if (particles_inside(i_particle)) then
                    do j_particle = i_particle + 1, num_particles
                        distance_ij = PBC_distance(Box_size, positions(:, i_particle), &
                                                             positions(:, j_particle))
                        if (distance_ij < cutoff) then
                            i_distribution = int(distance_ij/delta) + 1
                            distribution_step(i_distribution) = distribution_step(i_distribution) + 1._DP
                        end if
                    end do
                end if
            end do
            
            if (num_particles_step > 0) then
                distribution_step(:) = distribution_step(:) / real(num_particles_step, DP)
                distribution_function(:) = distribution_function(:) + distribution_step(:)
            end if
                
        end if

    end do
    call cpu_time(time_end)
    write(output_unit, *) "Finish !"

    deallocate(particles_inside)
    deallocate(positions)
    close(positions_unit)

    num_particles_inside = real(num_particles_sum, DP) / real(num_steps, DP)
    density_inside = num_particles_inside / Volume_inside

    open(newunit=distrib_unit, file=trim(name)//"_radial_distribution_function.out", &
         action="write")
    
        distribution_function(:) = 2._DP * distribution_function(:) / real(num_steps, DP)
    
        do i_distribution = 1, num_distribution
        
            distance_i = (real(i_distribution, DP) - 0.5_DP) * delta
            distance_i_minus = real(i_distribution - 1, DP) * delta
            distance_i_plus = real(i_distribution, DP) * delta
            
            distribution_function(i_distribution) = distribution_function(i_distribution) / &
                density_inside / (sphere_volume(distance_i_plus) - sphere_volume(distance_i_minus))
            write(distrib_unit, *) distance_i, distribution_function(i_distribution)
            
        end do
        
    close(distrib_unit)
    
    open(newunit=report_unit, file=trim(name)//"_radial_distribution_report.txt", &
         action="write")
         write(report_unit, *) "Volume inside:", Volume_inside
         write(report_unit, *) "Box lower bounds: ", Box_lower_bounds(:)
         write(report_unit, *) "Box upper bounds: ", Box_upper_bounds(:)
         write(report_unit, *) "Density inside: ", density_inside
         write(report_unit, *) "Cut off", cutoff
         write(report_unit, *) "Duration =", (time_end - time_start) / 60._DP, "min"
    close(report_unit)
    
    deallocate(distribution_function)
    deallocate(distribution_step)
    
    deallocate(domain_ratio)
    deallocate(Box_size)

end program radial_distribution
