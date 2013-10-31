!> \brief Calculate and print the distribution function

program distribution

use data_precisions, only : DP
use data_constants, only : PI
use data_box, only : Ndim, LsizeMi, Volume
use data_monteCarlo, only : Nstep
use data_distribution, only : snap, dist_dr
use module_physics, only : dist_PBC
use module_distribution, only : sphere_volume
!$ use omp_lib

implicit none

    character(len=5) :: name
    integer :: Ncol
    real(DP) :: rMin, rCut
    integer :: snap_factor
    real(DP) :: density
    integer, dimension(:), allocatable :: dist_sum
    integer :: snaps_unit, distrib_unit

    real(DP) :: rMax
    integer :: Ndist
    integer :: iStep
    integer :: iCol, jCol
    real(DP) :: r_ij
    integer :: iDist, iDistMin, iDistMax
    real(DP) :: r_iDist
    real(DP), dimension(:), allocatable :: dist_function
    real(DP), dimension(:, :), allocatable :: positions
    
    character(len=4096) :: file_name
    integer :: length, time_unit, file_stat

    real(DP) :: tIni, tFin
    !$ real(DP) :: tIni_para, tFin_para

    if (.not.snap) stop "Snap désactivé."
    
    call get_command_argument(1, file_name, length, file_stat)
    if (file_stat /= 0) stop "error get_command_argument"
    
    open(newunit=snaps_unit, recl=4096, file=file_name(1:length), status='old', action='read')
    
    read(snaps_unit, *) name, Ncol, rMin, rCut, snap_factor
    write(*, *) name, Ncol, rMin, rCut, snap_factor
    
    rMax = norm2(LsizeMi)
    Ndist = int(rMax/dist_dr)
    allocate(dist_sum(Ndist))
    allocate(dist_function(Ndist))
    allocate(positions(Ndim, Ncol))
    density = real(Ncol, DP) / Volume

    dist_sum(:) = 0   

    write(*, *) "Start !"
    call cpu_time(tIni)
    !$ tIni_para = omp_get_wtime()
    !$omp parallel do schedule(static) reduction(+:dist_sum) private(iCol, jCol, r_ij, iDist)
    do iStep = 1, Nstep/snap_factor

        ! Lecture :
        !$omp critical
        do iCol = 1, Ncol
            read(snaps_unit, *) positions(:, iCol)
        end do
        !$omp end critical

        ! Traitement
            
        do iCol = 1, Ncol
            do jCol = iCol + 1, Ncol

                r_ij = dist_PBC(positions(:, iCol), positions(:, jCol))      
                iDist =  int(r_ij/dist_dr)
                dist_sum(iDist) = dist_sum(iDist) + 1

            end do
        end do

    end do
    !$omp end parallel do
    !$ tFin_para = omp_get_wtime()
    call cpu_time(tFin)
    write(*, *) "Finish !"

    close(snaps_unit)

    ! Ecriture

    iDistMin = 0
    iDistMax = 0

    open(newunit=distrib_unit, file=name//"_dist_function.out", action="write")
    
        do iDist = 1, Ndist
        
            r_iDist = (real(iDist, DP) + 0.5_DP) * dist_dr
            dist_function(iDist) = 2._DP * real(dist_sum(iDist), DP) / real(Nstep/snap_factor, DP) / &
                real(Ncol, DP) / (sphere_volume(iDist+1)-sphere_volume(iDist)) / density
            write(distrib_unit, *) r_iDist, dist_function(iDist)
            
            if (r_iDist>=rMin .and. r_iDist<=rCut) then
                if (iDistMin == 0) then
                    iDistMin = iDist
                end if
                iDistMax = iDist
            end if
            
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
    deallocate(positions)

end program distribution
