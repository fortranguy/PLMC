!> \brief Distribution module

module mod_distrib

use data_constants
use data_distrib

implicit none

contains

    !> Calculate the volume of the sphere
    
    function sphereVol(iDist)
    
        integer, intent(in) :: iDist    
        real(DP) :: sphereVol
        
        sphereVol = 4._DP/3._DP * PI * ( real(iDist, DP)*deltaDist )**3
        
    end function sphereVol
    
end module mod_distrib

!> \brief Calculate and print the distribution function

program distribution

use data_constants
use data_distrib
use data_particles
use data_mc
use data_potentiel
use mod_distrib
use mod_physics
use class_interactingSpheres
!$ use omp_lib

implicit none

    real(DP), parameter :: densite = real(inter_Ncol, DP) / &
        (Lsize1*Lsize2*Lsize3)
    integer, dimension(:), allocatable :: distrib
    integer, parameter :: unitSnapShots = 10, unitDistrib = 11, &
        unitEnerg = 12

    integer :: iStep
    integer :: iCol, jCol
    real(DP) :: r_ij
    integer :: iDist, iDistMin, iDistMax
    real(DP) :: r
    real(DP) :: numerat, denomin
    real(DP), dimension(:), allocatable :: fct_dist
    real(DP) :: energSum	
    type(InteractingSpheres) :: inter
    real(DP), dimension(Dim, inter_Ncol) :: X

    !$ integer :: nb_taches
    real(DP) :: tIni, tFin
    !$ real(DP) :: tIni_para, tFin_para

    if (.not.snap) stop "Snap désactivé."

    call inter%construct()

    call initDistriParams()
    allocate(distrib(Ndist))
    allocate(fct_dist(Ndist))

    distrib(:) = 0

    open(unit=unitSnapShots, recl=4096, file="snap.shot", status='old', &
        action='read')

    call cpu_time(tIni)
    !$ tIni_para = omp_get_wtime()
    !$omp parallel private(X, iCol, jCol, r_ij, iDist)
    !$ nb_taches = omp_get_num_threads()
    !$omp do schedule(static, Nstep/nb_taches) reduction(+:distrib)
    do iStep = 1, Nstep

        ! Lecture :
        !$omp critical
        do iCol = 1, inter_Ncol
            read(unitSnapShots, *) X(:, iCol)
        end do
        !$omp end critical

        ! Traitement
            
        do iCol = 1, inter_Ncol
            do jCol = iCol + 1, inter_Ncol

                r_ij = dist(X(:, iCol), X(:, jCol))		        
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

    close(unitSnapShots)

    ! Ecriture

    iDistMin = 0
    iDistMax = 0

    open(unit=unitDistrib, file="fct_distrib.out", action="write")	
        do iDist = 1, Ndist
        
            r = (real(iDist, DP) + 0.5_DP) * deltaDist
            numerat = real(distrib(iDist), DP) / real(Nstep, DP)
            denomin = real(inter_Ncol, DP) * (sphereVol(iDist+1) - sphereVol(iDist))
            fct_dist(iDist) = 2._DP * numerat / denomin / densite
            write(unitDistrib, *) r, fct_dist(iDist)
            
            if (r>=inter_rMin .and. r<=inter_rCut) then
                if (iDistMin == 0) then
                    iDistMin = iDist
                end if
                iDistMax = iDist
            end if
            
        end do
    close(unitDistrib)

    ! Energie par particule

    energSum = 0._DP

    do iDist = iDistMin, iDistMax
        r = (real(iDist, DP) + 0.5_DP) * deltaDist
        energSum = energSum + inter%ePot_pair(r) * fct_dist(iDist) * 4._DP*PI*r**2	
    end do

    open(unit=unitEnerg, file="epp_dist.out", action="write")
        write(unitEnerg, *) "epp_dist =", &
            densite/2._DP * energSum * deltaDist
    close(unitEnerg)

    deallocate(fct_dist)
    deallocate(distrib)

end program distribution