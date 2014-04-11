!> \brief Calculate and write the distribution function

program distribution

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP
use data_box, only: Ndim, Box_size
use data_monte_carlo, only: Nstep
use data_distribution, only: snap, dist_dr
use module_physics_micro, only: sphere_volume, dist_PBC
!$ use omp_lib

implicit none

    character(len=5) :: name, name_bis
    integer :: num_particles, num_particles_bis
    integer :: snap_factor, snap_factor_bis
    real(DP) :: density
    integer, dimension(:), allocatable :: dist_sum
    integer :: positions_unit, orientations_unit, distrib_unit

    real(DP) :: rMax
    integer :: Ndist
    integer :: iStep
    integer :: i_particle, j_particle
    real(DP) :: r_ij
    integer :: iDist
    real(DP) :: r_iDist, r_iMinus, r_iPlus
    real(DP), dimension(:), allocatable :: dist_function
    real(DP), dimension(:, :), allocatable :: positions, orientations
    
    logical :: withOrientations
    
    character(len=4096) :: file
    integer :: length, time_unit, stat
    logical :: exist

    real(DP) :: tIni, tFin
    !$ real(DP) :: tIni_para, tFin_para

    if (.not.snap) stop "Snap désactivé."
    
    call get_command_argument(1, file, length, stat)
    if (stat /= 0) error stop "error get_command_argument"
    inquire(file=file(1:length), exist=exist)
    if (.not.exist) then
        write(error_unit, *) "missing file: ", file(1:length)
        error stop
    end if
    
    open(newunit=positions_unit, recl=4096, file=file(1:length), status='old', action='read')
    
    read(positions_unit, *) name, num_particles, snap_factor
    write(output_unit, *) name, num_particles, snap_factor
    
    rMax = norm2(Box_size/2._DP)
    Ndist = int(rMax/dist_dr)
    allocate(dist_sum(Ndist))
    allocate(dist_function(Ndist))
    allocate(positions(Ndim, num_particles))
    density = real(num_particles, DP) / product(Box_size)

    dist_sum(:) = 0
    
    withOrientations = (name == "dipol" .and. command_argument_count() == 2)
    
    if (withOrientations) then
    
        call get_command_argument(2, file, length, stat)
        if (stat /= 0) error stop "error get_command_argument"
        inquire(file=file(1:length), exist=exist)
        if (.not.exist) then
            write(error_unit, *) "missing file: ", file(1:length)
            error stop
        end if
        
        open(newunit=orientations_unit, recl=4096, file=file(1:length), status='old', &
        action='read')
        
        read(orientations_unit, *) name_bis, num_particles_bis, snap_factor_bis
        if (name_bis/=name .or. num_particles_bis/=num_particles .or. snap_factor_bis/=snap_factor) then
            write(error_unit, *) "Error: positions and orientations tags don't match."
            error stop
        end if
        
        allocate(orientations(Ndim, num_particles))
        
    end if

    write(output_unit, *) "Start !"
    call cpu_time(tIni)
    !$ tIni_para = omp_get_wtime()
    !$omp parallel do schedule(static) reduction(+:dist_sum) private(positions, i_particle, j_particle, &
    !$ r_ij, iDist)
    do iStep = 1, Nstep/snap_factor

        ! Read
        !$omp critical
        do i_particle = 1, num_particles
            read(positions_unit, *) positions(:, i_particle)
        end do
        
        if (withOrientations) then
            do i_particle = 1, num_particles
                read(orientations_unit, *) orientations(:, i_particle)
            end do
        end if
        !$omp end critical

        ! Fill
        do i_particle = 1, num_particles
            do j_particle = i_particle + 1, num_particles

                r_ij = dist_PBC(Box_size, positions(:, i_particle), positions(:, j_particle))
                iDist =  int(r_ij/dist_dr)
                dist_sum(iDist) = dist_sum(iDist) + 1

            end do
        end do

    end do
    !$omp end parallel do
    !$ tFin_para = omp_get_wtime()
    call cpu_time(tFin)
    write(output_unit, *) "Finish !"

    close(positions_unit)
    deallocate(positions)
    
    if (withOrientations) then
        close(orientations_unit)
        deallocate(orientations)
    end if

    open(newunit=distrib_unit, file=name//"_dist_function.out", action="write")
    
        do iDist = 1, Ndist
        
            r_iDist = (real(iDist, DP) + 0.5_DP) * dist_dr
            r_iMinus = real(iDist, DP) * dist_dr
            r_iPlus = real(iDist + 1, DP) * dist_dr
            
            dist_function(iDist) = 2._DP * real(dist_sum(iDist), DP) / real(Nstep/snap_factor, DP) / &
                real(num_particles, DP) / (sphere_volume(r_iPlus)-sphere_volume(r_iMinus)) / density
            write(distrib_unit, *) r_iDist, dist_function(iDist)
            
        end do
        
    close(distrib_unit)
    
    open(newunit=time_unit, file=name//"_dist_time.out")
        write(time_unit, *) "pseudo serial time", tFin - tIni
        !$ write(time_unit, *) "parallel time", tFin_para - tIni_para
        !$omp parallel
            !$omp master
                !$ write(time_unit, *) "number of threads =", omp_get_num_threads()
            !$omp end master
        !$omp end parallel
        !$ write(time_unit, *) "ratio =", (tFin-tIni)/(tFin_para-tIni_para)
            ! fake ?
    close(time_unit)
    
    deallocate(dist_function)
    deallocate(dist_sum)

end program distribution
