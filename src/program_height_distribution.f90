!> \brief Calculate and write the distribution function

!> 1. Density profile:
!> \f[
!>      \rho(z) = \frac{\langle N(z) \rangle \sigma^3}{S\delta z}
!> \f]

!> 2. Average orientation: different than II.40 ?
!> \f[
!>      Q_{zz}(z) = \left\langle \frac{\sum_{i=1}^N (3\mu_{i,z}^2/\mu_i^2 - 1)/2}
!>                                    {N(z)}\right\rangle
!> \f]

!> 3. Preferred orientation:
!> \f[
!>      Q = \frac{3}{2N} \sum_{i=1}^N |\vec{\mu}_i)(\vec{\mu}_i| - \frac{1}{2}I
!> \f]
!> We look for the eigen vector associated with the highest eigen value of Q.
!> \f[
!>      Q = \sum_{d=1}^3 |\vec{q}_d) q_d (\vec{q}_d|
!> \f]
!> \f[ \vec{d} := \vec{q_d} | q_d \text{max}\f]

program height_distribution

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit, error_unit
use data_precisions, only: real_zero
use data_box, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: data_filename, data_post_filename, report_filename, &
                       test_file_exists, test_data_found
use module_geometry, only: set_geometry, geometry
use module_linear_algebra, only: identity_matrix, eigen_symmetric
use module_physics_micro, only: PBC_distance
use module_arguments, only: arg_to_file

implicit none

    logical :: take_snapshot
    character(len=:), allocatable :: Box_geometry
    real(DP), dimension(:), allocatable :: Box_size
    real(DP) :: Box_height

    character(len=4096) :: positions_name, orientations_name
    integer :: num_particles, positions_num_particles, orientations_num_particles
    integer :: positions_snap_factor, orientations_snap_factor
    real(DP) :: density
    integer :: positions_unit, orientations_unit, distribution_unit
    
    integer :: num_thermalisation_steps
    integer :: num_equilibrium_steps, i_step, num_steps
    integer :: i_particle
    real(DP) :: distance_max, height_i, delta
    integer :: num_distribution, i_distribution
    real(DP), dimension(:, :), allocatable :: distribution_step
    real(DP), dimension(:, :), allocatable :: distribution_function
    real(DP), dimension(:, :), allocatable :: positions, orientations
    
    real(DP), dimension(num_dimensions) :: eigenvalues, preferred_orientation, &
                                           normed_orientation
    real(DP), dimension(num_dimensions, num_dimensions) :: orientation_matrix
    real(DP) :: orientation_z_sqr
    
    logical :: with_orientations
    
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

    call report_json%destroy()
    
    data_name = "Box.size"
    call data_json%get(data_name, Box_size, found)
    call test_data_found(data_name, found)
    
    if (geometry%bulk) then
        if (size(Box_size) /= num_dimensions) error stop "Box size dimension"
        Box_height = Box_size(3)
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

    call data_post_json%destroy()
    
    distance_max = Box_height
    num_distribution = int(distance_max/delta) + 1
    allocate(distribution_step(num_distribution, 3))
    allocate(distribution_function(num_distribution, 3))
    
    call arg_to_file(1, filename, length)
    open(newunit=positions_unit, recl=4096, file=filename(1:length), status='old', &
         action='read')
    
    read(positions_unit, *) positions_name, positions_num_particles, &
                            positions_snap_factor
    write(output_unit, *) trim(positions_name), positions_num_particles, &
                          positions_snap_factor
    
    num_particles = positions_num_particles
    allocate(positions(num_dimensions, num_particles))
    density = real(num_particles, DP) / product(Box_size)
    
    with_orientations = (command_argument_count() == 2)
    
    if (with_orientations) then
        
        call arg_to_file(2, filename, length)
        open(newunit=orientations_unit, recl=4096, file=filename(1:length), &
             status='old', action='read')
        
        read(orientations_unit, *) orientations_name, orientations_num_particles, &
                                   orientations_snap_factor
        if ((positions_name /= orientations_name) .or. &
            (positions_num_particles /= orientations_num_particles) .or. &
            (positions_snap_factor /= orientations_snap_factor)) then
            write(error_unit, *) "positions and orientations tags don't match."
            error stop
        end if
        
        allocate(orientations(num_dimensions, num_particles))
    
    end if

    write(output_unit, *) "Start !"
    call cpu_time(time_start)
    distribution_function(:, :) = 0._DP
    num_steps = 0
    do i_step = num_thermalisation_steps + 1, num_thermalisation_steps + num_equilibrium_steps
        if (modulo(i_step, positions_snap_factor) == 0) then
        
            num_steps = num_steps + 1
        
            do i_particle = 1, num_particles
                read(positions_unit, *) positions(:, i_particle)
            end do
            
            if (with_orientations) then
                do i_particle = 1, num_particles
                    read(orientations_unit, *) orientations(:, i_particle)
                end do
            end if

            ! Fill
            distribution_step(:, :) = 0
            do i_particle = 1, num_particles
                i_distribution = int(positions(3, i_particle)/delta) + 1
                distribution_step(i_distribution, 1) = distribution_step(i_distribution, 1) + 1._DP
            end do
            
            if (with_orientations) then
            
                orientation_matrix(:, :) = 0._DP
                do i_particle = 1, num_particles
                    orientation_matrix(:, :) = orientation_matrix(:, :) + &
                    matmul(reshape(orientations(:, i_particle), [num_dimensions, 1]), &
                           reshape(orientations(:, i_particle), [1, num_dimensions]))
                end do
                orientation_matrix(:, :) = 1.5_DP/real(num_particles, DP) * orientation_matrix(:, :) - &
                                           0.5_DP * identity_matrix(num_dimensions)
                call eigen_symmetric(orientation_matrix, eigenvalues)
                preferred_orientation(:) = orientation_matrix(:, num_dimensions) / &
                                           norm2(orientation_matrix(:, num_dimensions))
                                           
                do i_particle = 1, num_particles
                    i_distribution = int(positions(3, i_particle)/delta) + 1
                    orientation_z_sqr = orientations(3, i_particle)**2 / &
                                        dot_product(orientations(:, i_particle), &
                                                    orientations(:, i_particle))
                    distribution_step(i_distribution, 2) = distribution_step(i_distribution, 2) + &
                                                           1.5_DP*orientation_z_sqr - 0.5_DP
                    
                    normed_orientation(:) = orientations(:, i_particle) / &
                                            norm2(orientations(:, i_particle))
                    distribution_step(i_distribution, 3) = distribution_step(i_distribution, 3) + &
                                                           dot_product(normed_orientation, &
                                                                       preferred_orientation)
                end do
                where(distribution_step(:, 1) > real_zero)
                    distribution_step(:, 2) = distribution_step(:, 2) / distribution_step(:, 1)
                    distribution_step(:, 3) = abs(distribution_step(:, 3)) / distribution_step(:, 1)
                end where
            
            end if
            
            distribution_function(:, :) = distribution_function(:, :) + distribution_step(:, :)
        
        end if
    end do
    call cpu_time(time_end)
    write(output_unit, *) "Finish !"
    
    if (with_orientations) then
        close(orientations_unit)
        deallocate(orientations)
    end if

    close(positions_unit)
    deallocate(positions)
    
    distribution_function(:, :) = distribution_function(:, :) / real(num_steps, DP)
    
    open(newunit=report_unit, recl=4096, &
         file=trim(positions_name)//"_height_distribution_report.txt", action="write")
    write(report_unit, *) "Normalisation: ", sum(distribution_function(:, 1)) / real(num_particles, DP)
    write(report_unit, *) "With orientations: ", with_orientations
    write(report_unit, *) "Duration =", (time_end - time_start) / 60._DP, "min"
    close(report_unit)

    distribution_function(:, 1) = distribution_function(:, 1) / (Box_size(1)*Box_size(2)*delta)

    open(newunit=distribution_unit, recl=4096, &
         file=trim(positions_name)//"_height_distribution_function.out", action="write")
         
    if (with_orientations) then
        do i_distribution = 1, num_distribution
            height_i = (real(i_distribution, DP) - 0.5_DP) * delta
            write(distribution_unit, *) height_i, distribution_function(i_distribution, :)
        end do
    else
        do i_distribution = 1, num_distribution
            height_i = (real(i_distribution, DP) - 0.5_DP) * delta
            write(distribution_unit, *) height_i, distribution_function(i_distribution, 1)
        end do
    end if
        
    close(distribution_unit)
    
    deallocate(distribution_function)
    deallocate(distribution_step)
    
    deallocate(Box_size)

end program height_distribution
