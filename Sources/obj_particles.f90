!***********************************************************************
!* MODULE: Particles objects                                           *
!***********************************************************************

module obj_particles

use class_particle
use data_particles
use data_mc
use data_potentiel

implicit none

    type(Particle), protected :: sph
    
contains

    subroutine particle_init()
    
        ! Potential initializarion
        
        call ePotIni()
        
        ! Construction
                
        sph =   Particle(&        
                    radius = radius, &
                    rmin = rmin, &
                    Ncol = Ncol, &
                    dx = dx, &
                    rcut = rcut, &
                    pas = pas, &
                    iMin = iMin, &
                    Ntab = Ntab, &
                    epsilon = epsilon, &
                    alpha = alpha, &
                    Vtab = Vtab &
                )
        
    end subroutine particle_init

end module obj_particles
