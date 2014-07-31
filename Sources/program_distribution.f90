!> \brief Calculate and write the distribution function

program distribution

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit, error_unit
use data_box, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: test_data_found
use module_physics_micro, only: sphere_volume, PBC_distance
use module_arguments, only: arg_to_file

implicit none

    logical :: take_snapshot
    real(DP), dimension(:), allocatable :: Box_size
    

    character(len=4096) :: name
    integer :: num_particles
    integer :: snap_factor
    real(DP) :: density
    integer :: positions_unit, distrib_unit
    
    integer :: num_thermalisation_steps
    integer :: num_equilibrium_steps, i_step, num_steps
    integer :: i_particle, j_particle
    real(DP) :: distance_ij, distance_max, delta
    real(DP) :: distance_i_distribution, distance_i_minus, distance_i_plus
    integer :: num_distribution, i_distribution
    real(DP), dimension(:), allocatable :: distribution_step
    real(DP), dimension(:), allocatable :: distribution_function
    real(DP), dimension(:, :), allocatable :: positions
    
    type(json_file) :: json
    character(len=4096) :: data_name
    logical :: found
    character(len=4096) :: file_name
    integer :: length, time_unit

    real(DP) :: time_init, time_final
    
    call json_initialize()
    call json%load_file(filename = "data.json")
    
    data_name = "Distribution.take snapshot"
    call json%get(data_name, take_snapshot, found)
    call test_data_found(data_name, found)

    if (.not.take_snapshot) stop "No snap shots taken."
    
    data_name = "Box.size"
    call json%get(data_name, Box_size, found)
    call test_data_found(data_name, found)
    if (size(Box_size) /= num_dimensions) error stop "Box size dimension"
    
    data_name = "Monte Carlo.number of thermalisation steps"
    call json%get(data_name, num_thermalisation_steps, found)
    call test_data_found(data_name, found)
    
    data_name = "Monte Carlo.number of equilibrium steps"
    call json%get(data_name, num_equilibrium_steps, found)
    call test_data_found(data_name, found)
    
    data_name = "Distribution.delta"
    call json%get(data_name, delta, found)
    call test_data_found(data_name, found)
    
    call json%destroy()
    
    distance_max = norm2(Box_size / 2._DP)
    num_distribution = int(distance_max/delta)
    allocate(distribution_step(num_distribution))
    allocate(distribution_function(num_distribution))
    
    call arg_to_file(1, file_name, length)
    open(newunit=positions_unit, recl=4096, file=file_name(1:length), status='old', action='read')
    
    read(positions_unit, *) name, num_particles, snap_factor
    write(output_unit, *) trim(name), num_particles, snap_factor    
    
    allocate(positions(num_dimensions, num_particles))
    density = real(num_particles, DP) / product(Box_size)

    write(output_unit, *) "Start !"
    call cpu_time(time_init)
    distribution_function(:) = 0._DP
    num_steps = 0
    do i_step = num_thermalisation_steps + 1, num_thermalisation_steps + num_equilibrium_steps    
        if (modulo(i_step, snap_factor) == 0) then
        
            num_steps = num_steps + 1
        
            do i_particle = 1, num_particles
                read(positions_unit, *) positions(:, i_particle)
            end do

            ! Fill
            distribution_step(:) = 0
            do i_particle = 1, num_particles
                do j_particle = i_particle + 1, num_particles
                    distance_ij = PBC_distance(Box_size, positions(:, i_particle), &
                                                         positions(:, j_particle))
                    i_distribution =  int(distance_ij/delta)
                    distribution_step(i_distribution) = distribution_step(i_distribution) + 1._DP
                end do
            end do
            
            distribution_function(:) = distribution_function(:) + distribution_step(:)
        
        end if
    end do
    call cpu_time(time_final)
    write(output_unit, *) "Finish !"

    close(positions_unit)
    deallocate(positions)

    open(newunit=distrib_unit, file=trim(name)//"_distribution_function.out", action="write")
    
        distribution_function(:) = 2._DP * distribution_function(:) / real(num_steps, DP) / &
                                   real(num_particles, DP)
    
        do i_distribution = 1, num_distribution
        
            distance_i_distribution = (real(i_distribution, DP) + 0.5_DP) * delta
            distance_i_minus = real(i_distribution, DP) * delta
            distance_i_plus = real(i_distribution + 1, DP) * delta
            
            distribution_function(i_distribution) = distribution_function(i_distribution) / &
                density / (sphere_volume(distance_i_plus) - sphere_volume(distance_i_minus))
            write(distrib_unit, *) distance_i_distribution, distribution_function(i_distribution)
            
        end do
        
    close(distrib_unit)
    
    open(newunit=time_unit, file=trim(name)//"_distribution_duration.txt")
        write(time_unit, *) "duration", time_final - time_init
    close(time_unit)
    
    deallocate(distribution_function)
    deallocate(distribution_step)
    
    deallocate(Box_size)

end program distribution
