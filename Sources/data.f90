!> Base units
!>  u_length, u_moment, u_permittivity

!***************************************************************************************************
!> \brief Precisions data :
!> declaration of precisions parameters
!***************************************************************************************************
module data_precisions

use, intrinsic :: iso_fortran_env, only: REAL64

implicit none

    integer, parameter :: DP = REAL64 !< double precision
    real(DP), parameter :: real_zero = real(2**4, DP) * epsilon(1._DP)
    real(DP), parameter :: io_tiny = real(2**2, DP) * epsilon(1._DP)
    real(DP), parameter :: consist_tiny = real(2**13, DP) * epsilon(1._DP)

end module data_precisions
!***************************************************************************************************

!***************************************************************************************************
!> \brief Constants data :
!> declaration of constants
!***************************************************************************************************
module data_constants

use data_precisions, only: DP

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

use data_precisions, only: DP

implicit none
    
    integer, parameter :: Ndim = 3
    
    real(DP), parameter :: Box_size1 = 25._DP ! u_length
    real(DP), parameter :: Box_size2 = Box_size1 ! u_length
    real(DP), parameter :: Box_size3 = Box_size1 ! u_length
    real(DP), dimension(Ndim), parameter :: Box_size = [Box_size1, Box_size2, Box_size3] ! u_length

    integer, parameter :: Box_wave1 = 8 ! 1/u_length
    integer, parameter :: Box_wave2 = Box_wave1 ! 1/u_length
    integer, parameter :: Box_wave3 = Box_wave1 ! 1/u_length
    integer, dimension(Ndim), parameter :: Box_wave = [Box_wave1, Box_wave2, Box_wave3] ! 1/u_length
    
end module data_box
!***************************************************************************************************

!***************************************************************************************************
!> \brief Particles data :
!> declaration of the particles parameters
!***************************************************************************************************
module data_particles

use data_precisions, only: DP

implicit none

    integer, parameter :: dipol_num_particles = 281
    
    real(DP), parameter :: hard_diameter = 1._DP ! u_length
    integer, parameter :: hard_num_particles = 6750
    
    real(DP), parameter :: mix_delta = 0.2_DP ! u_length
    
end module data_particles
!***************************************************************************************************
    
!***************************************************************************************************
!> \brief Monte-Carlo data :
!> declaration of the Monte-Carlo parameters
!***************************************************************************************************
module data_monte_carlo

use data_precisions, only: DP
use data_box, only: Ndim

implicit none

    real(DP), parameter :: Temperature = 1._DP ! u_temperature
    integer, parameter :: decorrelFactor = 2**0
    integer, parameter :: switch_factor = 1
    integer, parameter :: Nthermal = 2**15/decorrelFactor
    integer, parameter :: Nadapt = 2**10
    integer, parameter :: Nstep = 2**16
    integer, parameter :: reset_iStep = 2**13/decorrelFactor
    
    ! move
    real(DP), dimension(Ndim), parameter :: dipol_move_delta = 0.3_DP ! u_length, adaptation
    real(DP), parameter :: dipol_move_rejectFix = 0.5_DP
    ! rotate
    real(DP), parameter :: dipol_rotate_delta = 30._DP ! u_moment, adaptation
    real(DP), parameter :: dipol_rotate_deltaMax = 75._DP
    real(DP), parameter :: dipol_rotate_rejectFix = 0.17_DP
    ! chemical potential
    integer, parameter :: dipol_Nwidom = 500
    
    real(DP), dimension(Ndim), parameter :: hard_move_delta = 0.5_DP ! u_length, adaptation
    real(DP), parameter :: hard_move_rejectFix = 0.5_DP
    integer, parameter :: hard_Nwidom = 500
    
end module data_monte_carlo
!***************************************************************************************************

!***************************************************************************************************
!> \brief Potential data :
!> declaration of the potential energy parameters

!> The dipolar spheres interaction uses the Ewald sums methods.

!> The mixing potential (mix) is composed of 3 parts :
!> hard sphere (HS) + Yukawa + cut.
!***************************************************************************************************
module data_potential

use data_precisions, only: DP

implicit none

    real(DP), parameter :: dipol_rMin_factor = 1._DP
    real(DP), parameter :: dipol_real_rCut_factor = 0.5_DP ! * Box_size(1)
    real(DP), parameter :: dipol_real_dr = 5.E-5_DP ! u_length
    real(DP), parameter :: dipol_alpha_factor = 7._DP ! / Box_size(1)

    real(DP), parameter :: hard_rMin_factor = 1._DP

    real(DP), parameter :: mix_rMin_factor = 1._dP
    real(DP), parameter :: mix_rCut = 1._DP ! u_length, adaptation
    real(DP), parameter :: mix_dr = 1._DP ! u_length
    real(DP), parameter :: mix_epsilon = 0._DP ! u_energy * u_length
    real(DP), parameter :: mix_alpha = 1._DP ! 1/u_length
    
    logical, parameter :: write_potential = .true.
        
end module data_potential
!***************************************************************************************************

!***************************************************************************************************
!> \brief Neighbours Cells data :
!> declaration of the grid/cell scheme parameters
!***************************************************************************************************
module data_neighbour_cells

use data_precisions, only: DP
use data_box, only: Ndim

implicit none

    integer, dimension(Ndim), parameter :: NnearCell_dim = 3 !< Number of nearest neighbour cells
                                                             !< in each direction
    integer, parameter :: NnearCell = NnearCell_dim(1) * NnearCell_dim(2) * NnearCell_dim(3)
                          !< Total number of nearest neighbour cells, including itself

end module data_neighbour_cells
!***************************************************************************************************

!***************************************************************************************************
!> \brief Histogram data :
!> declaration of the histogram function parameters
!***************************************************************************************************
module data_histogram

use data_precisions, only: DP

implicit none

    real(DP), parameter :: hist_dE = 1._DP

end module data_histogram
!***************************************************************************************************

!***************************************************************************************************
!> \brief Distribution data :
!> declaration of the distribution function parameters
!***************************************************************************************************
module data_distribution

use data_precisions, only: DP

implicit none

    logical, parameter :: snap = .true.
    real(DP), parameter :: dist_dr = 0.025_DP
    integer, parameter :: snap_ratio = 400

end module data_distribution
!***************************************************************************************************
