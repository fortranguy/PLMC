!***********************************************************************
!> \brief Constants data :
!> declaration of constants
!***********************************************************************
module data_constants

    ! double precision
    integer, parameter :: DP = selected_real_kind(15, 307)
    real(DP), parameter :: PI = acos(-1._DP)
        
end module data_constants
!***********************************************************************

!***********************************************************************
!> \brief Cell data :
!> declaration of the cell parameters
!***********************************************************************
module data_cell

use data_constants
    
    integer, parameter :: Dim = 3
    real(DP), parameter :: Lsize1 = 12._DP
    real(DP), parameter :: Lsize2 = Lsize1
    real(DP), parameter :: Lsize3 = Lsize1
    real(DP), dimension(Dim), parameter :: Lsize = &
        [Lsize1, Lsize2, Lsize3]
    real(DP), dimension(Dim), parameter :: LsizeMi = 0.5_DP * Lsize
    
end module data_cell
!***********************************************************************

!***********************************************************************
!> \brief Particles data :
!> declaration of the particles parameters
!***********************************************************************
module data_particles

use data_constants
use data_cell
    
    real(DP), parameter :: sph_radius = .5_DP
    real(DP), parameter :: sph_rmin = 1._DP
    integer, parameter :: sph_Ncol = 270
    
end module data_particles
!***********************************************************************
    
!***********************************************************************
!> \brief Monte-Carlo data :
!> declaration of the MC parameters
!***********************************************************************
module data_mc

use data_constants
use data_particles
use data_cell

implicit none

    real(DP), parameter :: Tstar = 1._DP
    integer, parameter :: Nstep = 2**10
    integer, parameter :: Ntherm = 2**8
    integer, parameter :: Nmove = 2**2 * sph_Ncol
    real(DP), dimension(Dim), parameter :: sph_dx = 2._DP

end module data_mc
!***********************************************************************

!***********************************************************************
!> \brief Potential data :
!> declaration of the potential energy parameters

!> The short potential is composed of 3 elements :
!> hard sphere (HS) + Yukawa + cut
!***********************************************************************
module data_potentiel

use data_constants
use data_particles

implicit none

    real(DP), parameter :: sph_rcut = 4._DP
    real(DP), parameter :: sph_pas = 5.E-5_DP
    integer, parameter :: sph_iMin = int(sph_rmin/sph_pas)
    integer, parameter :: sph_Ntab = int(sph_rcut/sph_pas)
    real(DP), parameter :: sph_epsilon = 1._DP
    real(DP), parameter :: sph_alpha = 5._DP
        
end module data_potentiel
!***********************************************************************

!***********************************************************************
!> \brief Neighbours data :
!> declaration of the grid/cell scheme parameters
!***********************************************************************
module data_neighbours

use data_cell

implicit none

    integer, dimension(Dim), parameter :: cell_neigh_coordMax = [3, 3, 3]
    integer, parameter :: cell_neighs_nb = 3**3 !< including itself

end module data_neighbours
!***********************************************************************

!***********************************************************************
!> \brief Distribution data :
!> declaration of the distribution function parameters
!***********************************************************************
module data_distrib

use data_constants
use data_cell

implicit none

	logical, parameter :: snap = .true.
	real(DP), parameter :: deltaDist = 0.01_DP
	real(DP), protected :: rMax
	integer, protected :: Ndist
	
contains
	
	subroutine initDistriParams()
	
	implicit none
	
		rMax = sqrt(dot_product(LsizeMi, LsizeMi))
		Ndist = int(rMax/deltaDist)
	
	end subroutine initDistriParams

end module data_distrib
!***********************************************************************