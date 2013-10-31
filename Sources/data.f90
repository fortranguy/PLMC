!> Base units
!>  u_length, u_moment, u_permittivity

!***************************************************************************************************
!> \brief Precisions data :
!> declaration of precisions paramters
!***************************************************************************************************
module data_precisions

use, intrinsic :: iso_fortran_env

implicit none

    integer, parameter :: DP = REAL64 !< double precision
    real(DP), parameter :: real_zero = real(2**3, DP) * epsilon(1._DP)
    real(DP), parameter :: io_tiny = real(2, DP) * epsilon(1._DP)
    real(DP), parameter :: consist_tiny = real(2**13, DP) * epsilon(1._DP)

end module data_precisions 
!***************************************************************************************************

!***************************************************************************************************
!> \brief Constants data :
!> declaration of constants
!***************************************************************************************************
module data_constants

use data_precisions, only : DP

implicit none
    
    real(DP), parameter :: PI = acos(-1._DP)
    real(DP), parameter :: sigma3d = 1._DP/sqrt(3._DP)    
        
end module data_constants
!***************************************************************************************************

!***************************************************************************************************
!> \brief Box data :
!> declaration of the simulation box parameters
!***************************************************************************************************
module data_box

use data_precisions, only : DP

implicit none
    
    integer, parameter :: Ndim = 3
    
    real(DP), parameter :: Lsize1 = 25._DP ! u_length
    real(DP), parameter :: Lsize2 = Lsize1 ! u_length
    real(DP), parameter :: Lsize3 = Lsize1 ! u_length 
    real(DP), dimension(Ndim), parameter :: Lsize = [Lsize1, Lsize2, Lsize3] ! u_length
    real(DP), dimension(Ndim), parameter :: LsizeMi = 0.5_DP * Lsize ! u_length
    real(DP), parameter :: Volume = Lsize(1) * Lsize(2) * Lsize(3) ! u_length**3

    integer, parameter :: Kmax1 = 8 ! 1/u_length
    integer, parameter :: Kmax2 = Kmax1 ! 1/u_length
    integer, parameter :: Kmax3 = Kmax1 ! 1/u_length
    integer, dimension(Ndim), parameter :: Kmax = [Kmax1, Kmax2, Kmax3] ! 1/u_length
    
    real(DP), parameter :: out_permittivity = 1._DP ! !< u_permittivity
    
end module data_box
!***************************************************************************************************

!***************************************************************************************************
!> \brief Particles data :
!> declaration of the particles parameters
!***************************************************************************************************
module data_particles

use data_precisions, only : DP

implicit none

    integer, parameter :: dipol_Ncol = 281
    
    real(DP), parameter :: hard_rMin = 1._DP ! u_length
    integer, parameter :: hard_Ncol = 6750
    
    real(DP), parameter :: mix_delta = 0.2_DP ! u_length
    
    ! out ---------------------------------------------------------------------
    real(DP), parameter :: inter_rMin = 1._DP ! u_length
    integer, parameter :: inter_Ncol = 0    
    ! -------------------------------------------------------------------------
    
end module data_particles
!***************************************************************************************************
    
!***************************************************************************************************
!> \brief Monte-Carlo data :
!> declaration of the Monte-Carlo parameters
!***************************************************************************************************
module data_monteCarlo

use data_precisions, only : DP
use data_box, only : Ndim

implicit none

    real(DP), parameter :: Temperature = 1._DP ! u_temperature
    integer, parameter :: decorrelFactor = 2**0
    integer, parameter :: Nthermal = 2**15/decorrelFactor
    integer, parameter :: Nadapt = 2**10
    integer, parameter :: Nstep = 2**16
    
    ! move
    real(DP), dimension(Ndim), parameter :: dipol_move_delta = 0.3_DP ! u_length, adaptation
    real(DP), parameter :: dipol_move_rejectFix = 0.5_DP
    ! rotate
    real(DP), parameter :: dipol_rotate_delta = 30._DP ! u_moment, adaptation
    real(DP), parameter :: dipol_rotate_deltaMax = 75._DP
    real(DP), parameter :: dipol_rotate_rejectFix = 0.17_DP
    ! chemical potential
    integer, parameter :: dipol_Nwidom = 500
    ! reinitialisations
    integer, parameter :: dipol_reInit_iStep = 2**13/decorrelFactor
    
    real(DP), dimension(Ndim), parameter :: hard_move_delta = 0.5_DP ! u_length, adaptation
    real(DP), parameter :: hard_move_rejectFix = 0.5_DP
    integer, parameter :: hard_Nwidom = 500
    
    ! out ---------------------------------------------------------------------
    real(DP), dimension(Ndim), parameter :: inter_move_delta = 1._DP ! u_length, adaptation
    real(DP), parameter :: inter_move_rejectFix = 0.5_DP
    integer, parameter :: inter_Nwidom = 0
    ! -------------------------------------------------------------------------

end module data_monteCarlo
!***************************************************************************************************

!***************************************************************************************************
!> \brief Potential data :
!> declaration of the potential energy parameters

!> The dipolar spheres interaction uses the Ewald sums methods.

!> The mixing potential (mix) is also composed of 3 parts :
!> hard sphere (HS) + Yukawa + cut.
!***************************************************************************************************
module data_potential

use data_precisions, only : DP
use data_box, only : Lsize1

implicit none

    real(DP), parameter :: dipol_rCut = Lsize1/2._DP ! u_length
    real(DP), parameter :: dipol_dr = 5.E-5_DP ! u_length
    real(DP), parameter :: dipol_alpha = 7._DP/Lsize1 ! 1/u_length
    
    real(DP), parameter :: mix_rCut = 1._DP ! u_length, adaptation
    real(DP), parameter :: mix_dr = 1._DP ! u_length
    real(DP), parameter :: mix_epsilon = 0._DP ! u_energy * u_length
    real(DP), parameter :: mix_alpha = 1._DP ! 1/u_length
    
    ! out ---------------------------------------------------------------------
    !> The interactive spheres (inter) potential is composed of 3 parts :
    !> hard sphere (HS) + Yukawa + cut.
    
    real(DP), parameter :: inter_rCut = 4._DP ! u_length
    real(DP), parameter :: inter_dr = 5.E-5_DP ! u_length
    real(DP), parameter :: inter_epsilon = 1._DP ! u_energy * u_length
    real(DP), parameter :: inter_alpha = 5._DP ! 1/u_length
    ! -------------------------------------------------------------------------
        
end module data_potential
!***************************************************************************************************

!***************************************************************************************************
!> \brief Neighbours Cells data :
!> declaration of the grid/cell scheme parameters
!***************************************************************************************************
module data_neighbourCells

use data_precisions, only : DP
use data_box, only : Ndim

implicit none

    integer, dimension(Ndim), parameter :: NnearCell_dim = 3 !< Number of nearest neighbour cells
                                                             !< in each direction
    integer, parameter :: NnearCell = NnearCell_dim(1) * NnearCell_dim(2) * NnearCell_dim(3)
                          !< Total number of nearest neighbour cells, including itself

end module data_neighbourCells
!***************************************************************************************************

!***************************************************************************************************
!> \brief Distribution data :
!> declaration of the distribution function parameters
!***************************************************************************************************
module data_distribution

use data_precisions, only : DP

implicit none

    logical, parameter :: snap = .true.
    real(DP), parameter :: deltaDist = 0.01_DP
    
    integer, parameter :: dipol_snap_factor = 1
    integer, parameter :: hard_snap_factor = 8
    
    ! out ---------------------------------------------------------------------
    integer, parameter :: inter_snap_factor = 1
    ! -------------------------------------------------------------------------

end module data_distribution
!***************************************************************************************************
