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

        ! Potential

        real(DP) :: rcut
        !real(DP) :: pas
        !real(DP) :: surpas
        !integer :: Ntab
        !integer :: iMin
        real(DP), dimension(:), allocatable :: Vtab
        !real(DP) :: epsilon
        !real(DP) :: alpha
        
    contains
    
        procedure :: rapport => Particle_rapport
        
    end type Particle
    
contains

    function Particle_rapport(this) result(rapport)
    
        class(Particle), intent(in) :: this
        real :: rapport        
        
        rapport = Lsize(1)/this%rayon
        
    end function Particle_rapport

end module class_particle

!***********************************************************************
!* MODULE: Particles objects                                           *
!***********************************************************************

module obj_particles

use class_particle
use data_particles

implicit none

    type(Particle), protected :: dip, sph
    
contains

    subroutine init()

        dip =   Particle(&
                    radius = ,&
                    rmin = ,&
                    Ncol = ,&
                    dx = ,&
                    rcut = ,&
                    Vtab = ,&
                )
                
        sph =   Particle(&
                    radius = ,&
                    rmin = ,&
                    Ncol = ,&
                    dx = ,&
                    rcut = ,&
                    Vtab = ,&
                )
        
    end subroutine init

end module obj_particles