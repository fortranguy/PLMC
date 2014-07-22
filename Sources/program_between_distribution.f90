!> \brief Calculate and print the distribution function between 2 types.

program between_distribution

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit, error_unit
use data_box, only: num_dimensions
use module_maths, only: lcm
use module_physics_micro, only: sphere_volume, PBC_distance
use module_arguments, only: arg_to_file
use module_read, only: jump_snap
!$ use omp_lib

implicit none
   
    logical :: take_snapshot

    character(len=5) :: type1_name, type2_name
    integer :: type1_num_particles, type2_num_particles
    integer :: type1_snap_factor, type2_snap_factor
    integer :: type1_positions_unit, type2_positions_unit
    integer :: snap_factors_lcm
    
    integer :: report_unit, distrib_unit
    integer :: num_distribution, i_distribution
    integer :: num_steps, i_step
    integer :: type1_i_particle, type2_i_particle
    real(DP) :: distance_12
    real(DP) :: distance_i, distance_minus, distance_plus
    real(DP), parameter :: distance_max = 3._DP
    real(DP), dimension(:), allocatable :: distribution_step
    real(DP), dimension(:), allocatable :: distribution_function
    real(DP), dimension(:, :), allocatable :: type1_positions, type2_positions
    
    character(len=4096) :: file
    integer :: length
    
    real(DP) :: initial_time, final_time
    
    if (.not.take_snapshot) stop "No snap shots taken."
    
    Volume_inside = product(Lsize(1:2)) * (z_max_ratio - z_min_ratio) * (Height-1._DP)
    
    num_distribution = int(distance_max / dist_dr)
    allocate(distribution_step(num_distribution))
    allocate(distribution_function(num_distribution))
    
    call arg_to_file(1, file, length)
    open(newunit=type1_positions_unit, recl=4096, file=file(1:length), status='old', action='read')    
    read(type1_positions_unit, *) type1_name, type1_num_particles, type1_snap_factor
    write(output_unit, *) "type 1: ", type1_name, type1_num_particles, type1_snap_factor
    allocate(type1_positions(Ndim, type1_num_particles))
    allocate(type1_paricles_inside(type1_num_particles))
    
    call arg_to_file(2, file, length)
    open(newunit=type2_positions_unit, recl=4096, file=file(1:length), status='old', action='read')    
    read(type2_positions_unit, *) type2_name, type2_num_particles, type2_snap_factor
    write(output_unit, *) "type 1: ", type2_name, type2_num_particles, type2_snap_factor
    allocate(type2_positions(Ndim, type2_num_particles))
    allocate(type2_paricles_inside(type2_num_particles))
    
    snap_factors_lcm = lcm(type1_snap_factor, type2_snap_factor)
    
    type1_num_particles_sum = 0
    type2_num_particles_sum = 0
    
    write(output_unit, *) "Start !"
    call cpu_time(initial_time)
    do i_step = 1, num_steps / snap_factors_lcm
        
        call jump_snap(snap_factors_lcm, type1_snap_factor, type1_num_particles, type1_positions_unit)
        do type1_i_particle = 1, type1_num_particles
            read(type1_positions_unit, *) type1_positions(:, type1_i_particle)
        end do

        call jump_snap(snap_factors_lcm, type2_snap_factor, type2_num_particles, type2_positions_unit)
        do type2_i_particle = 1, type2_num_particles
            read(type2_positions_unit, *) type2_positions(:, type2_i_particle)
        end do
        
        call test_particle_inside(type1_num_particles, type1_positions, type1_paricles_inside, type1_num_particles_step)
        type1_num_particles_sum = type1_num_particles_sum + type1_num_particles_step
        call test_particle_inside(type2_num_particles, type2_positions, type2_paricles_inside, type2_num_particles_step)
        type2_num_particles_sum = type2_num_particles_sum + type2_num_particles_step
        
        distribution_step(:) = 0._DP
        
        do type1_i_particle = 1, type1_num_particles
            if (type1_paricles_inside(type1_i_particle))  then
            
                do type2_i_particle = 1, type2_num_particles                
                    distance_12 = PBC_distance(type1_positions(:, type1_i_particle), type2_positions(:, type2_i_particle))
                    if (distance_12 <= distance_max) then
                        i_distribution = int(distance_12 / dist_dr)
                        distribution_step(i_distribution) = distribution_step(i_distribution) + 1._DP
                    end if
                end do
                
            end if
        end do
        
        distribution_step(:) = distribution_step(:) / real(type1_num_particles_step, DP)
        distribution_function(:) = distribution_function(:) + distribution_step(:)
    
    end do
    call cpu_time(final_time)
    write(output_unit, *) "Finish !"
    
    deallocate(type2_paricles_inside)
    deallocate(type2_positions) 
    close(type2_positions_unit)   
    
    deallocate(type1_paricles_inside)
    deallocate(type1_positions) 
    close(type1_positions_unit)
    
    open(newunit=report_unit, file=type1_name//"-"//type2_name//"_mix_distribution_report.txt", &
         action="write")
    
    write(report_unit, *) "Volume_inside =", Volume_inside
    write(report_unit, *) "z-domain = ", z_min_ratio * Height, z_max_ratio * Height
    
    type1_num_particles_inside = real(type1_num_particles_sum, DP) / real(num_steps/snap_factors_lcm, DP)
    type1_density_inside = type1_num_particles_inside / Volume_inside
    write(report_unit, *) type1_name, " inside density: ", type1_density_inside
    type2_num_particles_inside = real(type2_num_particles_sum, DP) / real(num_steps/snap_factors_lcm, DP)
    type2_density_inside = type2_num_particles_inside / Volume_inside
    write(report_unit, *) type2_name, " inside density: ", type2_density_inside
    write(report_unit, *) "Duration =", (final_time - initial_time) / 60._DP, "min"
    
    close(report_unit)
    
    open(newunit=distrib_unit, file=type1_name//"-"//type2_name//"_mix_distribution.out", &
         action="write")
    
        distribution_function(:) = distribution_function(:) / real(num_steps/snap_factors_lcm, DP)
    
        do i_distribution = 1, num_distribution
        
            distance_i = (real(i_distribution, DP) + 0.5_DP) * dist_dr
            distance_minus = real(i_distribution, DP) * dist_dr
            distance_plus = real(i_distribution + 1, DP) * dist_dr
            
            distribution_function(i_distribution) = distribution_function(i_distribution) / type2_density_inside / &
                                   (sphere_volume(distance_plus) - sphere_volume(distance_minus))
            write(distrib_unit, *) distance_i, distribution_function(i_distribution)
            
        end do
        
    close(distrib_unit)
    
    deallocate(distribution_function)
    deallocate(distribution_step)

end program between_distribution
