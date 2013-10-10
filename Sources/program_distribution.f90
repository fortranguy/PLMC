!> \brief Distribution module

module module_distrib

use data_precisions, only : DP
use data_constants, only : PI
use data_distribution, only : deltaDist

implicit none

contains

    !> Calculate the volume of the sphere
    
    pure function sphereVol(iDist)
    
        integer, intent(in) :: iDist    
        real(DP) :: sphereVol
        
        sphereVol = 4._DP/3._DP * PI * (real(iDist, DP)*deltaDist)**3
        
    end function sphereVol
    
end module module_distrib

!> \brief Calculate and print the distribution function

program distribution

use data_precisions, only : DP
use data_constants, only : PI
use data_box, only : LsizeMi, Volume
use data_particles
use data_monteCarlo
use data_potential
use data_distribution
use module_physics
use module_distrib
use class_hardSpheres
use class_interactingSpheres
use class_dipolarSpheres
use class_mixingPotential
use module_algorithms
!$ use omp_lib

implicit none

    real(DP), parameter :: densite = real(hard_Ncol, DP) / Volume
    integer, dimension(:), allocatable :: distrib
    integer, parameter :: snaps_unit = 10, distrib_unit = 11, energ_unit = 12

    real(DP) :: rMax
    integer :: Ndist
    integer :: iStep
    integer :: iCol, jCol
    real(DP) :: r_ij
    integer :: iDist, iDistMin, iDistMax
    real(DP) :: r
    real(DP) :: numerat, denomin
    real(DP), dimension(:), allocatable :: fct_dist
    real(DP) :: energSum
    type(HardSpheres) :: hard
    type(MixingPotential) :: mix
    real(DP), dimension(Ndim, hard_Ncol) :: X

    !$ integer :: nb_taches
    real(DP) :: tIni, tFin
    !$ real(DP) :: tIni_para, tFin_para

    if (.not.snap) stop "Snap désactivé."

    call mix%construct()
    call hard%construct()

    rMax = norm2(LsizeMi)
    Ndist = int(rMax/deltaDist)
    allocate(distrib(Ndist))
    allocate(fct_dist(Ndist))

    distrib(:) = 0

    open(unit=snaps_unit, recl=4096, file="snap.shot", status='old', action='read')

    call cpu_time(tIni)
    !$ tIni_para = omp_get_wtime()
    !$omp parallel private(X, iCol, jCol, r_ij, iDist)
    !$ nb_taches = omp_get_num_threads()
    !$omp do schedule(static, Nstep/nb_taches) reduction(+:distrib)
    do iStep = 1, Nstep

        ! Lecture :
        !$omp critical
        do iCol = 1, hard_Ncol
            read(snaps_unit, *) X(:, iCol)
        end do
        !$omp end critical

        ! Traitement
            
        do iCol = 1, hard_Ncol
            do jCol = iCol + 1, hard_Ncol

                r_ij = dist_PBC(X(:, iCol), X(:, jCol))      
                iDist =  int(r_ij/deltaDist)
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
        
            r = (real(iDist, DP) + 0.5_DP) * deltaDist
            numerat = real(distrib(iDist), DP) / real(Nstep, DP)
            denomin = real(hard_Ncol, DP) * (sphereVol(iDist+1)-sphereVol(iDist))
            fct_dist(iDist) = 2._DP * numerat / denomin / densite
            write(distrib_unit, *) r, fct_dist(iDist)
            
            if (r>=hard%get_rMin() .and. r<=hard%get_rCut()) then
                if (iDistMin == 0) then
                    iDistMin = iDist
                end if
                iDistMax = iDist
            end if
            
        end do
    close(distrib_unit)

    ! Energie par particule

    energSum = 0._DP

    do iDist = iDistMin, iDistMax
        r = (real(iDist, DP) + 0.5_DP) * deltaDist
        energSum = energSum + hard%Epot_pair(r) * fct_dist(iDist) * 4._DP*PI*r**2
    end do

    open(unit=energ_unit, file="epp_dist.out", action="write")
        write(energ_unit, *) "epp_dist =", densite/2._DP * energSum * deltaDist
    close(energ_unit)

    deallocate(fct_dist)
    deallocate(distrib)

    call hard%destroy()
    call mix%destroy()

end program distribution
