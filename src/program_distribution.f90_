!> \brief Calculate and print the distribution function

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

program distribution

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP, real_zero
use data_box, only: Ndim, Lsize, Height
use data_monteCarlo, only: Nstep
use data_distribution, only: snap, dist_dz
use mkl95_lapack, only: SYEV
!$ use omp_lib

implicit none

    character(len=5) :: name, nameBis
    integer :: Ncol, NcolBis
    integer :: snap_factor, snap_factorBis
    integer :: positions_unit, orientations_unit, distrib_unit

    integer :: i_dim
    integer :: Ndist
    integer :: iStep
    integer :: iCol
    integer :: iDist
    real(DP) :: z_iDist
    real(DP) :: mCol_z_sqr
    real(DP), dimension(:, :), allocatable :: dist_step
    real(DP), dimension(:, :), allocatable :: dist_function
    real(DP), dimension(:, :), allocatable :: positions, orientations
    real(DP), dimension(Ndim) :: eigenvalues, preferred_orientation, normed_orientation
    real(DP), dimension(Ndim, Ndim) :: orientation_matrix, identity_matrix
    
    logical :: withOrientations
    
    character(len=4096) :: filename
    integer :: length, time_unit, stat
    logical :: exist

    real(DP) :: tIni, tFin
    !$ real(DP) :: tIni_para, tFin_para

    if (.not.snap) stop "Snap désactivé."
    
    call get_command_argument(1, filename, length, stat)
    if (stat /= 0) error stop "error get_command_argument"
    inquire(file=filename(1:length), exist=exist)
    if (.not.exist) then
        write(error_unit, *) "missing file: ", filename(1:length)
        error stop
    end if
    
    open(newunit=positions_unit, recl=4096, file=filename(1:length), status='old', action='read')
    
    read(positions_unit, *) name, Ncol, snap_factor
    write(output_unit, *) name, Ncol, snap_factor

    Ndist = int(Height/dist_dz)
    allocate(dist_step(Ndist, 3))
    allocate(dist_function(Ndist, 3))
    allocate(positions(Ndim, Ncol))
    
    withOrientations = (name == "dipol" .and. command_argument_count() == 2)
    
    if (withOrientations) then
    
        call get_command_argument(2, filename, length, stat)
        if (stat /= 0) error stop "error get_command_argument"
        inquire(file=filename(1:length), exist=exist)
        if (.not.exist) then
            write(error_unit, *) "missing file: ", filename(1:length)
            error stop
        end if
        
        open(newunit=orientations_unit, recl=4096, file=filename(1:length), status='old', &
        action='read')
        
        read(orientations_unit, *) nameBis, NcolBis, snap_factorBis
        if (nameBis/=name .or. NcolBis/=Ncol .or. snap_factorBis/=snap_factor) then
            write(error_unit, *) "Error: positions and orientations tags don't match."
            error stop
        end if
        
        allocate(orientations(Ndim, Ncol))      
        
    end if

    dist_function(:, :) = 0._DP
    identity_matrix(:, :) = 0._DP
    forall (i_dim = 1:Ndim) identity_matrix(i_dim, i_dim) = 1._DP

    write(output_unit, *) "Start !"
    call cpu_time(tIni)
    !$ tIni_para = omp_get_wtime()
    !$omp parallel do schedule(static) reduction(+:dist_step, dist_function) &
    !$                                 private(positions, orientations, iCol, iDist, mCol_z_sqr)
    do iStep = 1, Nstep/snap_factor
    
        dist_step(:, :) = 0._DP

        ! Read
        !$omp critical
        do iCol = 1, Ncol
            read(positions_unit, *) positions(:, iCol)
        end do
        
        if (withOrientations) then
            do iCol = 1, Ncol
                read(orientations_unit, *) orientations(:, iCol)
            end do
        end if
        !$omp end critical

        ! Fill
        do iCol = 1, Ncol
            iDist =  int(positions(3, iCol)/dist_dz)
            dist_step(iDist, 1) = dist_step(iDist, 1) + 1._DP
        end do

        if (withOrientations) then
                
            orientation_matrix(:, :) = 0._DP
            do iCol = 1, Ncol
                orientation_matrix(:, :) = orientation_matrix(:, :) + &
                                           matmul(reshape(orientations(:, iCol), [Ndim, 1]), &
                                                  reshape(orientations(:, iCol), [1, Ndim]))
            end do
            orientation_matrix(:, :) = 1.5_DP/real(Ncol, DP) * orientation_matrix(:, :) - &
                                       0.5_DP * identity_matrix(:, :)
            call syev(orientation_matrix, eigenvalues, 'V')
            preferred_orientation(:) = orientation_matrix(:, Ndim) / norm2(orientation_matrix(:, Ndim))
            
            do iCol = 1, Ncol
                iDist =  int(positions(3, iCol)/dist_dz)
                
                mCol_z_sqr = orientations(3, iCol)**2 / dot_product(orientations(:, iCol), &
                                                                    orientations(:, iCol))
                dist_step(iDist, 2) = dist_step(iDist, 2) + 1.5_DP*mCol_z_sqr - 0.5_DP
                
                normed_orientation(:) = orientations(:, iCol) / norm2(orientations(:, iCol))
                dist_step(iDist, 3) = dist_step(iDist, 3) + dot_product(normed_orientation, &
                                                                        preferred_orientation)
            end do                                    
            where(dist_step(:, 1) > real_zero)
                dist_step(:, 2) = dist_step(:, 2) / dist_step(:, 1)
                dist_step(:, 3) = abs(dist_step(:, 3)) / dist_step(:, 1)
            end where
            
        end if

        dist_function(:, :) = dist_function(:, :) + dist_step(:, :)

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
    
    dist_function(:, :) = dist_function(:, :) / real(Nstep/snap_factor, DP)

    write(output_unit, *) "checking normalisation"
    write(output_unit, *) sum(dist_function(:, 1)) / real(Ncol, DP)

    dist_function(:, 1) = dist_function(:, 1) / (Lsize(1)*Lsize(2)*dist_dz)
    
    open(newunit=distrib_unit, recl=4096, file=name//"_dist_function.out", action="write")
    
    if (withOrientations) then
        do iDist = 1, Ndist
            z_iDist = (real(iDist, DP) + 0.5_DP) * dist_dz
            write(distrib_unit, *) z_iDist, dist_function(iDist, :)
        end do
    else
        do iDist = 1, Ndist
            z_iDist = (real(iDist, DP) + 0.5_DP) * dist_dz
            write(distrib_unit, *) z_iDist, dist_function(iDist, 1)
        end do
    end if
    
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
    deallocate(dist_step)

end program distribution
