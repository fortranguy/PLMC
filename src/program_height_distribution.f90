!> \brief Calculate and write the distribution function

program height_distribution

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit, error_unit
use data_box, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: test_data_found
use module_geometry, only: set_geometry
use module_physics_micro, only: PBC_distance
use module_arguments, only: arg_to_file

implicit none

    logical :: take_snapshot
    real(DP), dimension(:), allocatable :: Box_size    
    real(DP) :: Box_height

    character(len=4096) :: positions_name, orientations_name
    integer :: positions_num_particles, orientations_num_particles
    integer :: positions_snap_factor, orientations_snap_factor
    real(DP) :: density
    integer :: positions_unit, orientations_unit, distrib_unit
    
    integer :: num_thermalisation_steps
    integer :: num_equilibrium_steps, i_step, num_steps
    integer :: i_particle
    real(DP) :: distance_max, height_i, delta
    integer :: num_distribution, i_distribution
    real(DP), dimension(:), allocatable :: distribution_step
    real(DP), dimension(:), allocatable :: distribution_function
    real(DP), dimension(:, :), allocatable :: positions, orientations
    
    logical :: with_orientations
    
    type(json_file) :: data_json, report_json
    character(len=4096) :: data_name
    logical :: found
    character(len=4096) :: filename
    integer :: length
    integer :: report_unit
    real(DP) :: time_start, time_end

    character(len=:), allocatable :: geometry
    
    call json_initialize()
    call data_json%load_file(filename = "data.json")
    
    data_name = "Distribution.take snapshot"
    call data_json%get(data_name, take_snapshot, found)
    call test_data_found(data_name, found)

    if (.not.take_snapshot) stop "No snap shots taken."

    call report_json%load_file(filename = "report.json")
    
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
    
    if (geometry == "bulk") then
        Box_height = Box_size(3)
    else if (geometry == "slab") then
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
    
    data_name = "Distribution.delta"
    call data_json%get(data_name, delta, found)
    call test_data_found(data_name, found)
    
    call data_json%destroy()
    
    distance_max = Box_height
    num_distribution = int(distance_max/delta)
    allocate(distribution_step(num_distribution))
    allocate(distribution_function(num_distribution))
    
    call arg_to_file(1, filename, length)
    open(newunit=positions_unit, recl=4096, file=filename(1:length), status='old', &
         action='read')
    
    read(positions_unit, *) positions_name, positions_num_particles, &
                            positions_snap_factor
    write(output_unit, *) trim(positions_name), positions_num_particles, &
                          positions_snap_factor    
    
    allocate(positions(num_dimensions, positions_num_particles))
    density = real(positions_num_particles, DP) / product(Box_size)
    
    with_orientations = (command_argument_count() == 2)
    
    if (with_orientations) then
        
        call arg_to_file(2, filename, length)
        open(newunit=orientations_unit, recl=4096, file=filename(1:length), &
             status='old', action='read')
        
        read(positions_unit, *) orientations_name, orientations_num_particles, &
                                orientations_snap_factor
        if ((positions_name /= orientations_name) .or. &
            (positions_num_particles /= orientations_num_particles) .or. &
            (positions_snap_factor /= orientations_snap_factor)) then
            write(error_unit, *) "positions and orientations tags don't match."
            error stop
        end if
        
        allocate(orientations(num_dimensions, orientations_num_particles))
    
    end if

    write(output_unit, *) "Start !"
    call cpu_time(time_start)
    distribution_function(:) = 0._DP
    num_steps = 0
    do i_step = num_thermalisation_steps + 1, num_thermalisation_steps + num_equilibrium_steps    
        if (modulo(i_step, positions_snap_factor) == 0) then
        
            num_steps = num_steps + 1
        
            do i_particle = 1, positions_num_particles
                read(positions_unit, *) positions(:, i_particle)
            end do

            ! Fill
            distribution_step(:) = 0
            do i_particle = 1, positions_num_particles
                do j_particle = i_particle + 1, positions_num_particles
                    distance_ij = PBC_distance(Box_size, positions(:, i_particle), &
                                                         positions(:, j_particle))
                    i_distribution =  int(distance_ij/delta)
                    distribution_step(i_distribution) = distribution_step(i_distribution) + 1._DP
                end do
            end do
            
            distribution_function(:) = distribution_function(:) + distribution_step(:)
        
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

    open(newunit=distrib_unit, file=trim(positions_name)//"_height_distribution_function.out", &
         action="write")
    
        distribution_function(:) = 2._DP * distribution_function(:) / real(num_steps, DP) / &
                                   real(positions_num_particles, DP)
    
        do i_distribution = 1, num_distribution
        
            distance_i_distribution = (real(i_distribution, DP) + 0.5_DP) * delta
            distance_i_minus = real(i_distribution, DP) * delta
            distance_i_plus = real(i_distribution + 1, DP) * delta
            
            distribution_function(i_distribution) = distribution_function(i_distribution) / &
                density / (sphere_volume(distance_i_plus) - sphere_volume(distance_i_minus))
            write(distrib_unit, *) distance_i_distribution, distribution_function(i_distribution)
            
        end do
        
    close(distrib_unit)
    
    open(newunit=report_unit, file=trim(positions_name)//"_height_distribution_report.txt", &
         action="write")
        write(report_unit, *) "Duration =", (time_end - time_start) / 60._DP, "min"
    close(report_unit)
    
    deallocate(distribution_function)
    deallocate(distribution_step)
    
    deallocate(Box_size)

end program height_distribution
