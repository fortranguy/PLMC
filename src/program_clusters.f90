program cluster

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit, error_unit
use data_box, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: data_filename, data_post_filename, report_filename, &
                       test_file_exists, test_data_found
use module_geometry, only: set_geometry, geometry
use module_physics_micro, only: PBC_vector, dipolar_pair_energy
use module_clusters, only: pairs_to_clusters
use module_arguments, only: arg_to_file

implicit none

    logical :: take_snapshot
    character(len=:), allocatable :: Box_geometry
    real(DP), dimension(:), allocatable :: Box_size
    integer :: num_steps

    character(len=4096) :: name, name_bis
    integer :: num_particles, num_particles_bis
    integer :: i_particle, j_particle
    integer :: snap_factor, snap_factor_bis
    integer :: positions_unit, orientations_unit
    
    character(len=4096) :: file
    integer :: length
    
    integer :: i_step
    real(DP), dimension(:, :), allocatable :: all_positions, all_orientations
    real(DP) :: cutoff
    
    real(DP) :: pair_energy
    real(DP), dimension(num_dimensions) :: vector_ij
    logical :: ij_linked

    integer :: num_sizes, i_size
    integer :: cluster_size
    integer :: num_clusters, i_cluster
    integer :: num_pairs
    integer, dimension(:), allocatable :: clusters_sizes
    integer, dimension(:, :), allocatable :: all_pairs
    real(DP), dimension(:), allocatable :: clusters_sizes_distribution
    real(DP) :: delta
    integer :: clusters_distribution_unit

    logical :: with_orientations
    
    type(json_file) :: data_json, data_post_json, report_json
    character(len=4096) :: data_name
    logical :: found
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
    else if (geometry%slab) then
        if (size(Box_size) /= 2) error stop "Box size dimension"
    end if
    
    data_name = "Monte Carlo.number of equilibrium steps"
    call data_json%get(data_name, num_steps, found)
    call test_data_found(data_name, found)

    call data_json%destroy()

    call test_file_exists(data_post_filename)
    call data_post_json%load_file(filename = data_post_filename)

    data_name = "Clusters.cut off"
    call data_post_json%get(data_name, cutoff, found)
    call test_data_found(data_name, found)
    
    data_name = "Clusters.number of sizes"
    call data_post_json%get(data_name, num_sizes, found)
    call test_data_found(data_name, found)    
    
    call data_post_json%destroy()

    call arg_to_file(1, file, length)
    open(newunit=positions_unit, recl=4096, file=file(1:length), status='old', action='read')
    read(positions_unit, *) name, num_particles, snap_factor
    write(output_unit, *) trim(name), num_particles, snap_factor
    allocate(all_positions(num_dimensions, num_particles))
    num_pairs = num_particles * (num_particles - 1) / 2
    delta = real(num_sizes, DP) / real(num_particles, DP)
    allocate(clusters_sizes(num_pairs))
    allocate(all_pairs(2, num_pairs))
    allocate(clusters_sizes_distribution(0:num_sizes))

    with_orientations = (command_argument_count() == 2)

    if (with_orientations) then
    
        call arg_to_file(2, file, length)
        open(newunit=orientations_unit, recl=4096, file=file(1:length), status='old', &
            action='read')
        read(orientations_unit, *) name_bis, num_particles_bis, snap_factor_bis
        if (name_bis/=name .or. num_particles_bis/=num_particles .or. snap_factor_bis/=snap_factor) then
            write(error_unit, *) "Error: positions and orientations tags don't match."
            error stop
        end if
        allocate(all_orientations(num_dimensions, num_particles))

    end if
    
    clusters_sizes_distribution(:) = 0
    
    clusters_sizes(:) = 0
    
    write(output_unit, *) "Start !"
    call cpu_time(time_start)
    do i_step = 1, num_steps/snap_factor
    
        do i_particle = 1, num_particles
            read(positions_unit, *) all_positions(:, i_particle)
        end do
        if (with_orientations) then
            do i_particle = 1, num_particles
                read(orientations_unit, *) all_orientations(:, i_particle)
            end do
        end if

        num_pairs = 0
        do i_particle = 1, num_particles
            do j_particle = i_particle + 1, num_particles

                ij_linked = .false.
                vector_ij(:) = PBC_vector(Box_size, &
                                          all_positions(:, i_particle), all_positions(:, j_particle))
                if (norm2(vector_ij) < cutoff) then
                    if (with_orientations) then
                        pair_energy = dipolar_pair_energy(all_orientations(:, i_particle), &
                                                        all_orientations(:, j_particle), &
                                                        vector_ij)
                        if (pair_energy < 0._DP) then
                            ij_linked = .true.
                        end if
                    else
                        ij_linked = .true.
                    end if
                    if (ij_linked) then
                        num_pairs = num_pairs + 1
                        all_pairs(:, num_pairs) = [i_particle, j_particle]
                    end if
                end if
                                                                
            end do
        end do
        
        call pairs_to_clusters(num_particles, all_pairs, num_pairs, clusters_sizes, num_clusters)
        
        do i_cluster = 1, num_clusters
            cluster_size = clusters_sizes(i_cluster)
            i_size = floor(real(cluster_size, DP) * delta)
            clusters_sizes_distribution(i_size) = clusters_sizes_distribution(i_size) + 1._DP
        end do
    
    end do
    call cpu_time(time_end)
    write(output_unit, *) "Finish !"
    
    clusters_sizes_distribution(0:num_sizes) = clusters_sizes_distribution(0:num_sizes) / &
                                               real(num_steps/snap_factor, DP)
    clusters_sizes_distribution(0:num_sizes) = clusters_sizes_distribution(0:num_sizes) / &
sum(clusters_sizes_distribution(0:num_sizes))
    
    open(newunit=clusters_distribution_unit, &
         file=trim(name)//"_clusters_sizes_distribution_histogram.out", action="write")
    do i_size = 1, num_sizes
        write(clusters_distribution_unit, *) real(i_size, DP) / real(num_sizes, DP), &
                                             clusters_sizes_distribution(i_size) / delta
                                             
    end do
    close(clusters_distribution_unit)

    open(newunit=report_unit, file=trim(name)//"_clusters_sizes_distribution_report.txt", &
         action="write")
        write(report_unit, *) "With orientations: ", with_orientations
        write(report_unit, *) "Duration =", (time_end - time_start) / 60._DP, "min"
    close(report_unit)

    if (with_orientations) deallocate(all_orientations)
    deallocate(clusters_sizes_distribution)
    deallocate(all_pairs)
    deallocate(clusters_sizes)
    deallocate(all_positions)
    
    if (with_orientations) close(orientations_unit)
    close(positions_unit)
    
    deallocate(Box_size)

end program cluster
