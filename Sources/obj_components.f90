!***********************************************************************
!* MODULE: Particles components                                        *
!***********************************************************************

module obj_components

use class_component
use data_particles
use data_mc
use data_potentiel

implicit none

    type(Component), protected :: sph
    
contains

    subroutine component_init()
    
        ! Component initializarion
        
        call ePotIni()
        
        ! Construction
        
        sph =   Component(&        
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
        
    end subroutine component_init

end module obj_components
