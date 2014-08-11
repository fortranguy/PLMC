program cluster

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit, error_unit
use data_box, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: test_data_found
use module_physics_micro, only: PBC_vector, dipolar_pair_energy
use module_clusters, only: pairs_to_clusters
use module_arguments, only: arg_to_file

implicit none

    logical :: take_snapshot
    real(DP), dimension(:), allocatable :: Box_size
    integer :: num_steps

    character(len=5) :: name, name_bis
    integer :: num_particles, num_particles_bis
    integer :: i_particle, j_particle
    integer :: snap_factor, snap_factorBis
    integer :: positions_unit, orientations_unit
    
    character(len=4096) :: file
    integer :: length 
    
    integer :: i_step
    real(DP), dimension(:, :), allocatable :: all_positions, all_orientations
    real(DP), parameter :: range_cut = 1.3_DP
    
    real(DP) :: pair_energy
    real(DP), dimension(num_dimensions) :: vector_ij

    integer :: cluster_size, cluster_size_max
    integer :: num_clusters, i_cluster
    integer :: num_pairs
    integer, dimension(:), allocatable :: clusters_sizes
    integer, dimension(:, :), allocatable :: all_pairs
    real(DP), dimension(:), allocatable :: clusters_sizes_distribution
    integer :: clusters_distribution_unit
    
    type(json_file) :: data_json
    character(len=4096) :: data_name
    logical :: found
    
    call json_initialize()
    call data_json%load_file(filename = "data.json")
    
    data_name = "Distribution.take snapshot"
    call data_json%get(data_name, take_snapshot, found)
    call test_data_found(data_name, found)
    
    if (.not.take_snapshot) stop "No snap shots taken."
    
    data_name = "Box.size"
    call data_json%get(data_name, Box_size, found)
    call test_data_found(data_name, found)
    if (size(Box_size) /= num_dimensions) error stop "Box size dimension"
    
    data_name = "Monte Carlo.number of equilibrium steps"
    call data_json%get(data_name, num_steps, found)
    call test_data_found(data_name, found)
    
    call data_json%destroy()

    call arg_to_file(1, file, length)    
    open(newunit=positions_unit, recl=4096, file=file(1:length), status='old', action='read')    
    read(positions_unit, *) name, num_particles, snap_factor
    write(output_unit, *) name, num_particles, snap_factor
    allocate(all_positions(num_dimensions, num_particles))
    num_pairs = num_particles * (num_particles - 1) / 2
    allocate(clusters_sizes(num_pairs))
    allocate(all_pairs(2, num_pairs))
    allocate(clusters_sizes_distribution(num_pairs))
    
    call arg_to_file(2, file, length) 
    open(newunit=orientations_unit, recl=4096, file=file(1:length), status='old', &
        action='read')    
    read(orientations_unit, *) name_bis, num_particles_bis, snap_factorBis
    if (name_bis/=name .or. num_particles_bis/=num_particles .or. snap_factorBis/=snap_factor) then
        write(error_unit, *) "Error: positions and orientations tags don't match."
        error stop
    end if
    allocate(all_orientations(num_dimensions, num_particles))
    
    clusters_sizes_distribution(:) = 0
    cluster_size_max = 0
    
    clusters_sizes(:) = 0
    
    write(output_unit, *) "Start !"
    do i_step = 1, num_steps/snap_factor
    
        do i_particle = 1, num_particles
            read(positions_unit, *) all_positions(:, i_particle)
        end do        
        do i_particle = 1, num_particles
            read(orientations_unit, *) all_orientations(:, i_particle)
        end do

        num_pairs = 0
        do i_particle = 1, num_particles
            do j_particle = i_particle + 1, num_particles
                
                vector_ij(:) = PBC_vector(Box_size, &
                                          all_positions(:, i_particle), all_positions(:, j_particle))
                if (norm2(vector_ij) < range_cut) then
                    pair_energy = dipolar_pair_energy(all_orientations(:, i_particle), &
                                                      all_orientations(:, j_particle), &
                                                      vector_ij)
                    if (pair_energy < 0._DP) then
                        num_pairs = num_pairs + 1
                        all_pairs(:, num_pairs) = [i_particle, j_particle]                        
                    end if
                end if
                                                                
            end do
        end do
        
        call pairs_to_clusters(num_particles, all_pairs, num_pairs, clusters_sizes, num_clusters)
        
        do i_cluster = 1, num_clusters
            cluster_size = clusters_sizes(i_cluster)
            clusters_sizes_distribution(cluster_size) = clusters_sizes_distribution(cluster_size) + 1._DP
            if (cluster_size_max < cluster_size) cluster_size_max = cluster_size
        end do
    
    end do
    write(output_unit, *) "Finish !"
    
    clusters_sizes_distribution(1:cluster_size_max) = clusters_sizes_distribution(1:cluster_size_max) / &
                                                     real(num_steps/snap_factor, DP)
    
    open(newunit=clusters_distribution_unit, file="clusters_sizes_distribution.out", action="write")
    do cluster_size = 1, cluster_size_max
        write(clusters_distribution_unit, *) cluster_size, clusters_sizes_distribution(cluster_size)
    end do
    close(clusters_distribution_unit)
    
    deallocate(all_orientations)
    deallocate(clusters_sizes_distribution)
    deallocate(all_pairs)
    deallocate(clusters_sizes)
    deallocate(all_positions)
    
    close(orientations_unit)
    close(positions_unit)
    
    deallocate(Box_size)

end program cluster
