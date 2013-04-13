!> \brief Monte Carlo simulation in canonical ensemble for a mixture

program mc_canonical

use data_constants
use data_mc
use data_distrib
use mod_tools
use class_interacting
use class_observables
use class_units

implicit none

! Beginning

    ! Initialisation
    
    integer :: iStep, iMove
    type(Interacting) :: inter_sph
    type(Observables) :: inter_obs
    type(Units) :: inter_io
    
    real(DP) :: tIni, tFin
        
    inter_sph = inter_constructor()
    call inter_obs%init()
    call inter_io%open()
    
    write(*, *) "Monte-Carlo - Canonical : Volume =", product(Lsize)    
    
    call inter_sph%report(inter_io%report)
    call init_random_seed(inter_io%report)
    
    ! Initial condition
    
    call initialCondition(inter_io%report, inter_sph%X)
    call inter_sph%snapShot(inter_io%snapIni)
    
    call inter_sph%overlapTest()
    inter_obs%ePot_total = inter_sph%ePot_total()
    call inter_sph%cols_to_cells()
    
! Middle --------------------------------------------------

    write(*, *) "Beginning of cycles"
    
    call cpu_time(tIni)
    do iStep = 1, Ntherm + Nstep
    
        do iMove = 1, Nmove
            call inter_sph%mcMove(inter_obs%ePot_total, inter_obs%Nrejects)
        end do
        
        call inter_sph%widom(inter_obs%activExInv)
        
        if (iStep <= Ntherm) then
        
            call inter_sph%adapt_dx(iStep, inter_obs%rejectsRateSum, inter_io%report)
            write(inter_io%dx, *) iStep, inter_sph%get_dx(), &
                inter_obs%rejectsRateSum/real(iStep, DP)
            write(inter_io%obsTherm, *) iStep, inter_obs%ePot_total, &
                inter_obs%activExInv
        
        else
            
            call inter_obs%add()
            write(inter_io%obs, *) iStep, inter_obs%ePot_total, &
                inter_obs%activExInv
            
            if (snap) then
                call inter_sph%snapShot(inter_io%snapShots)
            end if
            
        end if
        
        call inter_obs%addReject()
    
    end do
    call cpu_time(tFin)

    write(*, *) "End of cycles"

! End -----------------------------------------------------

    call inter_sph%overlapTest()
    
    call inter_sph%consistTest(inter_obs%ePot_total, inter_io%report)    
    call inter_obs%results(tFin-tIni, inter_io%report)    
    call inter_sph%snapShot(inter_io%snapFin)
    
    call inter_sph%destructor()
    call inter_io%close()
    
end program mc_canonical