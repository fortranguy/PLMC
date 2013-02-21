!***********************************************************************
!* MODULE: Particle class                                              *
!***********************************************************************

module class_particle

use data_cell

implicit none

private

    type, public :: Particle

        ! Particles

        real(DP) :: radius
        real(DP) :: rmin
        integer ::  Ncol

        ! Monte-Carlo
        
        real(DP), dimension(Dim) :: dx

        ! Potential domain

        real(DP) :: rcut
        real(DP) :: pas
        integer :: iMin
        integer :: Ntab
        
        ! Potential shape
        
        real(DP) :: epsilon
        real(DP) :: alpha
        real(DP), dimension(:), allocatable :: Vtab
        
    contains
    
        procedure :: ePot
        
    end type Particle
    
contains

    function ePot(this, r)
        
        class(Particle), intent(in) :: this
        real, intent(in) :: r
        
        integer :: i
        real :: r_i, ePot
       
        if (r < this%rcut) then
       
            i = int(r/this%pas)
            r_i = real(i, DP)*this%pas
            ePot = this%Vtab(i) + (r-r_i)/this%pas * &
                (this%Vtab(i+1)-this%Vtab(i))
           
        else
       
            ePot = 0.
           
        end if
        
    end function ePot

end module class_particle

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