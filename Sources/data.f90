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
    
    integer, parameter :: num_dimensions = 3
    
end module data_box
!***************************************************************************************************

!***************************************************************************************************
!> \brief Neighbours Cells data :
!> declaration of the grid/cell scheme parameters
!***************************************************************************************************
module data_neighbour_cells

use data_precisions, only: DP
use data_box, only: num_dimensions

implicit none

    integer, dimension(num_dimensions), parameter :: num_near_cells_dim = 3 !< Number of nearest
                                                                            !< neighbour cells
                                                                            !< in each direction
    integer, parameter :: num_near_cells = num_near_cells_dim(1) * num_near_cells_dim(2) * &
                                           num_near_cells_dim(3) !< Total number of nearest neighbour
                                                                 !< cells, including itself

end module data_neighbour_cells
!***************************************************************************************************
