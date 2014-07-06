!> \brief Calculate and write the distribution function

program distribution

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP
use data_box, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: test_data_found
use module_physics_micro, only: sphere_volume, PBC_distance
use module_arguments, only: arg_to_file
!$ use omp_lib

implicit none

    logical :: take_snapshot
    real(DP), dimension(:), allocatable :: Box_size
    integer :: num_steps

    character(len=4096) :: name, name_bis
    integer :: num_particles, num_particles_bis
    integer :: snap_factor, snap_factor_bis
    real(DP) :: density
    integer, dimension(:), allocatable :: distribution_sum
    integer :: positions_unit, orientations_unit, distrib_unit
    
    real(DP) :: max_distance, delta
    integer :: i_step
    integer :: i_particle, j_particle
    real(DP) :: distance_ij
    real(DP) :: distance_i_distribution, distance_i_minus, distance_i_plus
    integer :: num_distribution, i_distribution
    real(DP), dimension(:), allocatable :: distribution_function
    real(DP), dimension(:, :), allocatable :: positions, orientations
    
    logical :: with_orientations
    
    type(json_file) :: json
    character(len=4096) :: data_name
    logical :: found
    character(len=4096) :: file_name
    integer :: length, time_unit

    real(DP) :: time_init, time_final
    !$ real(DP) :: time_init_para, time_final_para
    
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
    
    data_name = "Monte Carlo.number of equilibrium steps"
    call json%get(data_name, num_steps, found)
    call test_data_found(data_name, found)
    
    data_name = "Distribution.delta"
    call json%get(data_name, delta, found)
    call test_data_found(data_name, found)
    
    call json%destroy()
    
    call arg_to_file(1, file_name, length)
    open(newunit=positions_unit, recl=4096, file=file_name(1:length), status='old', action='read')
    
    read(positions_unit, *) name, num_particles, snap_factor
    write(output_unit, *) trim(name), num_particles, snap_factor
    
    max_distance = norm2(Box_size/2._DP)
    num_distribution = int(max_distance/delta)
    allocate(distribution_sum(num_distribution))
    allocate(distribution_function(num_distribution))
    allocate(positions(num_dimensions, num_particles))
    density = real(num_particles, DP) / product(Box_size)

    distribution_sum(:) = 0
    
    with_orientations = (command_argument_count() == 2)
    
    if (with_orientations) then
    
        call arg_to_file(2, file_name, length)
        open(newunit=orientations_unit, recl=4096, file=file_name(1:length), status='old', &
        action='read')
        
        read(orientations_unit, *) name_bis, num_particles_bis, snap_factor_bis
        if (name_bis/=name .or. num_particles_bis/=num_particles .or. snap_factor_bis/=snap_factor) then
            write(error_unit, *) "Error: positions and orientations tags don't match."
            error stop
        end if
        
        allocate(orientations(num_dimensions, num_particles))
        
    end if

    write(output_unit, *) "Start !"
    call cpu_time(time_init)
    !$ time_init_para = omp_get_wtime()
    !$omp parallel do schedule(static) reduction(+:distribution_sum) &
    !$ private(positions, i_particle, j_particle, distance_ij, i_distribution)
    do i_step = 1, num_steps/snap_factor

        ! Read
        !$omp critical
        do i_particle = 1, num_particles
            read(positions_unit, *) positions(:, i_particle)
        end do
        
        if (with_orientations) then
            do i_particle = 1, num_particles
                read(orientations_unit, *) orientations(:, i_particle)
            end do
        end if
        !$omp end critical

        ! Fill
        do i_particle = 1, num_particles
            do j_particle = i_particle + 1, num_particles

                distance_ij = PBC_distance(Box_size, positions(:, i_particle), positions(:, j_particle))
                i_distribution =  int(distance_ij/delta)
                distribution_sum(i_distribution) = distribution_sum(i_distribution) + 1

            end do
        end do

    end do
    !$omp end parallel do
    !$ time_final_para = omp_get_wtime()
    call cpu_time(time_final)
    write(output_unit, *) "Finish !"

    close(positions_unit)
    deallocate(positions)
    
    if (with_orientations) then
        close(orientations_unit)
        deallocate(orientations)
    end if

    open(newunit=distrib_unit, file=trim(name)//"_distribution_function.out", action="write")
    
        do i_distribution = 1, num_distribution
        
            distance_i_distribution = (real(i_distribution, DP) + 0.5_DP) * delta
            distance_i_minus = real(i_distribution, DP) * delta
            distance_i_plus = real(i_distribution + 1, DP) * delta
            
            distribution_function(i_distribution) = &
                2._DP * real(distribution_sum(i_distribution), DP) / real(num_steps / &
                snap_factor, DP) / real(num_particles, DP) / &
                (sphere_volume(distance_i_plus) - sphere_volume(distance_i_minus)) / density
            write(distrib_unit, *) distance_i_distribution, distribution_function(i_distribution)
            
        end do
        
    close(distrib_unit)
    
    open(newunit=time_unit, file=trim(name)//"_distribution_time.out")
        write(time_unit, *) "pseudo serial time", time_final - time_init
        !$ write(time_unit, *) "parallel time", time_final_para - time_init_para
        !$omp parallel
            !$omp master
                !$ write(time_unit, *) "number of threads =", omp_get_num_threads()
            !$omp end master
        !$omp end parallel
        !$ write(time_unit, *) "ratio =", (time_final-time_init)/(time_final_para-time_init_para)
            ! fake ?
    close(time_unit)
    
    deallocate(distribution_function)
    deallocate(distribution_sum)
    
    deallocate(Box_size)

end program distribution
