!> \brief Calculate and print the distribution function

program distribution

use data_precisions, only : DP
use data_constants, only : PI
use data_box, only : Ndim, LsizeMi, Volume
use data_monteCarlo, only : Nstep
use data_distribution, only : snap, dist_dr
use module_physics, only : dist_PBC
use module_distribution, only : sphereVol
!$ use omp_lib

implicit none

    character(len=5) :: name
    integer :: Ncol
    integer :: rMin, rCut
    real(DP) :: densite
    integer, dimension(:), allocatable :: distrib
    integer :: snaps_unit, distrib_unit = 11

    real(DP) :: rMax
    integer :: Ndist
    integer :: iStep
    integer :: iCol, jCol
    real(DP) :: r_ij
    integer :: iDist, iDistMin, iDistMax
    real(DP) :: r
    real(DP) :: numerat, denomin
    real(DP), dimension(:), allocatable :: fct_dist
    real(DP), dimension(:, :), allocatable :: positions
    
    character(len=20) :: file_name
    integer :: length, file_stat

    !$ integer :: nb_taches
    real(DP) :: tIni, tFin
    !$ real(DP) :: tIni_para, tFin_para

    if (.not.snap) stop "Snap désactivé."
    
    call get_command_argument(1, file_name, length, file_stat)
    if (file_stat /= 0) stop "error get_command_argument"
    open(newunit=snaps_unit, recl=4096, file=file_name(1:length), status='old', action='read')
    
    read(snaps_unit, *) name, Ncol, rMin, rCut
    write(*, *) name, Ncol, rMin, rCut
    stop
    
    rMax = norm2(LsizeMi)
    Ndist = int(rMax/dist_dr)
    allocate(distrib(Ndist))
    allocate(fct_dist(Ndist))
    allocate(positions(Ndim, Ncol))
    densite = real(Ncol, DP) / Volume

    distrib(:) = 0

    

    call cpu_time(tIni)
    !$ tIni_para = omp_get_wtime()
    !$omp parallel private(positions, iCol, jCol, r_ij, iDist)
    !$ nb_taches = omp_get_num_threads()
    !$omp do schedule(static, Nstep/nb_taches) reduction(+:distrib)
    do iStep = 1, Nstep

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
                distrib(iDist) = distrib(iDist) + 1

            end do
        end do

    end do
    !$omp end do nowait
    !$omp end parallel
    call cpu_time(tFin)
    !$ tFin_para = omp_get_wtime()

    open(unit=100, file="dist_duree.out")
        write(100, *) "DuréeSérie_pseudo", tFin - tIni
        !$ write(100, *) "DuréeParallèle", tFin_para - tIni_para
        !$ write(100, *) "nb_taches =", nb_taches
        !$ write(100, *) "Rapport =", (tFin-tIni)/(tFin_para-tIni_para)
            ! trompeur ?
    close(100)

    close(snaps_unit)

    ! Ecriture

    iDistMin = 0
    iDistMax = 0

    open(unit=distrib_unit, file="fct_distrib.out", action="write")
        do iDist = 1, Ndist
        
            r = (real(iDist, DP) + 0.5_DP) * dist_dr
            numerat = real(distrib(iDist), DP) / real(Nstep, DP)
            denomin = real(Ncol, DP) * (sphereVol(iDist+1)-sphereVol(iDist))
            fct_dist(iDist) = 2._DP * numerat / denomin / densite
            write(distrib_unit, *) r, fct_dist(iDist)
            
            if (r>=rMin .and. r<=rCut) then
                if (iDistMin == 0) then
                    iDistMin = iDist
                end if
                iDistMax = iDist
            end if
            
        end do
    close(distrib_unit)
    
    deallocate(fct_dist)
    deallocate(distrib)
    deallocate(positions)

end program distribution
