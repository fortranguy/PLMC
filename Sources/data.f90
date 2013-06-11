!***************************************************************************************************
!> \brief Constants data :
!> declaration of constants
!***************************************************************************************************
module data_constants

implicit none

    ! double precision
    integer, parameter :: DP = selected_real_kind(15, 307)
    real(DP), parameter :: PI = acos(-1._DP)
    real(DP), parameter :: sigma3d = 1._DP/sqrt(3._DP)
    real(DP), parameter :: io_tiny = 1.E-10_DP
    real(DP), parameter :: consist_tiny = 1.E-10_DP
        
end module data_constants
!***************************************************************************************************

!***************************************************************************************************
!> \brief Cell data :
!> declaration of the cell parameters
!***************************************************************************************************
module data_cell

use data_constants

implicit none
    
    integer, parameter :: Dim = 3
    
    real(DP), parameter :: Lsize1 = 25._DP
    real(DP), parameter :: Lsize2 = Lsize1
    real(DP), parameter :: Lsize3 = Lsize1
    real(DP), dimension(Dim), parameter :: Lsize = [Lsize1, Lsize2, Lsize3]
    real(DP), dimension(Dim), parameter :: LsizeMi = 0.5_DP * Lsize
    real(DP), parameter :: Volume = Lsize1 * Lsize2 * Lsize3

    integer, parameter :: Kmax1 = 8
    integer, parameter :: Kmax2 = Kmax1
    integer, parameter :: Kmax3 = Kmax1
    integer, dimension(Dim), parameter :: Kmax = [Kmax1, Kmax2, Kmax3]
    
end module data_cell
!***************************************************************************************************

!***************************************************************************************************
!> \brief Particles data :
!> declaration of the particles parameters
!***************************************************************************************************
module data_particles

use data_constants
use data_cell

implicit none

    real(DP), parameter :: dipol_radius = .5_DP
    real(DP), parameter :: dipol_rMin = 2._DP * dipol_radius
    integer, parameter :: dipol_Ncol = 281

    real(DP), parameter :: inter_radius = .5_DP
    real(DP), parameter :: inter_rMin = 2._DP * inter_radius
    integer, parameter :: inter_Ncol = 0
    
    real(DP), parameter :: hard_radius = .5_DP
    real(DP), parameter :: hard_rMin = 2._DP * hard_radius
    integer, parameter :: hard_Ncol = 6750
    
    real(DP), parameter :: mix_delta = 0.2_DP
    real(DP), parameter :: mix_rMin = dipol_radius + hard_radius + mix_delta
    
    integer, parameter :: Ncol = dipol_Ncol + hard_Ncol
    
end module data_particles
!***************************************************************************************************
    
!***************************************************************************************************
!> \brief Monte-Carlo data :
!> declaration of the Monte-Carlo parameters
!***************************************************************************************************
module data_mc

use data_constants
use data_particles
use data_cell

implicit none

    real(DP), parameter :: Tstar = 1._DP
    integer, parameter :: Nstep = 2**16
    integer, parameter :: decorrelFactor = 2**0
    integer, parameter :: Ntherm = 2**15/decorrelFactor
    integer, parameter :: Nmove = decorrelFactor * Ncol
    integer, parameter :: Nrotate = decorrelFactor * dipol_Ncol
    
    integer, parameter :: dipol_structure_iStep = 2**13/decorrelFactor
    ! move
    real(DP), dimension(Dim), parameter :: dipol_deltaX = 0.3_DP
    real(DP), parameter :: dipol_rejectFix = 0.5_DP
    integer, parameter :: dipol_Nadapt = Ntherm/50
    ! rotate
    real(DP), parameter :: dipol_deltaM = 30._DP
    real(DP), parameter :: dipol_deltaMmax = 75._DP
    real(DP), parameter :: dipol_rejectRotFix = 0.17_DP
    integer, parameter :: dipol_NadaptRot = Ntherm/50
    ! chemical potential
    integer, parameter :: dipol_Nwidom = 500
    
    real(DP), dimension(Dim), parameter :: inter_deltaX = 1._DP
    real(DP), parameter :: inter_rejectFix = 0.5_DP
    integer, parameter :: inter_Nadapt = Ntherm/8
    integer, parameter :: inter_Nwidom = inter_Ncol
    
    real(DP), dimension(Dim), parameter :: hard_deltaX = 0.5_DP
    real(DP), parameter :: hard_rejectFix = 0.5_DP
    integer, parameter :: hard_Nadapt = Ntherm/50
    integer, parameter :: hard_Nwidom = 500

end module data_mc
!***************************************************************************************************

!***************************************************************************************************
!> \brief Potential data :
!> declaration of the potential energy parameters

!> The interactive spheres (inter) potential is composed of 3 parts :
!> hard sphere (HS) + Yukawa + cut

!> The mixing potential (mix) is also composed of 3 parts :
!> hard sphere (HS) + Yukawa + cut
!***************************************************************************************************
module data_potentiel

use data_constants
use data_particles

implicit none

    real(DP), parameter :: dipol_rCut = Lsize1/2._DP*sqrt(3._DP)
    real(DP), parameter :: dipol_dr = 5.E-5_DP
    real(DP), parameter :: dipol_alpha = 7._DP/Lsize1

    real(DP), parameter :: inter_rCut = 4._DP
    real(DP), parameter :: inter_dr = 5.E-5_DP
    real(DP), parameter :: inter_epsilon = 1._DP
    real(DP), parameter :: inter_alpha = 5._DP
    
    real(DP), parameter :: hard_rCut = hard_rMin
    
    real(DP), parameter :: mix_rCut = 1.25_DP
    real(DP), parameter :: mix_dr = mix_rCut/2._DP
    real(DP), parameter :: mix_epsilon = 0._DP
    real(DP), parameter :: mix_alpha = 40._DP
        
end module data_potentiel
!***************************************************************************************************

!***************************************************************************************************
!> \brief Neighbours data :
!> declaration of the grid/cell scheme parameters
!***************************************************************************************************
module data_neighbours

use data_constants
use data_cell
use data_potentiel

implicit none

    integer, dimension(Dim), parameter :: cell_neigh_coordMax = 3
    integer, parameter :: cell_neighs_nb = 3**3 !< including itself

    real(DP), dimension(Dim), parameter :: dipol_cell_Lsize = dipol_rMin
    real(DP), dimension(Dim), parameter :: inter_cell_Lsize = inter_rCut
    real(DP), dimension(Dim), parameter :: hard_cell_Lsize = hard_rCut
    
    real(DP), dimension(Dim), parameter :: mix_cell_Lsize = mix_rCut

end module data_neighbours
!***************************************************************************************************

!***************************************************************************************************
!> \brief Distribution data :
!> declaration of the distribution function parameters
!***************************************************************************************************
module data_distrib

use data_constants
use data_cell

implicit none

    logical, parameter :: snap = .true.
    real(DP), parameter :: deltaDist = 0.01_DP
    real(DP), protected :: rMax
    integer, protected :: Ndist
    
    integer, parameter :: dipol_snap_factor = 1
    integer, parameter :: inter_snap_factor = 1
    integer, parameter :: hard_snap_factor = 8

contains

    subroutine initDistriParams()

        rMax = norm2(LsizeMi)
        Ndist = int(rMax/deltaDist)

    end subroutine initDistriParams

end module data_distrib
!***************************************************************************************************
