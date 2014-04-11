!> \brief Definition of derives types

module module_types

use data_precisions, only: DP
use data_box, only: Ndim

implicit none

    ! System Box
    
    type, public :: Box_dimensions
        real(DP), dimension(Ndim) :: size
        integer, dimension(Ndim) :: wave
    end type Box_dimensions

    ! For neighbour cells

    type, public :: Node
        integer :: number
        type(Node), pointer :: next => null()
    end type Node

    type, public :: LinkedList
        type(Node), pointer :: particle => null()
    end type LinkedList
    
    ! For algorithms
    
    type, public :: particle_index
        integer :: number
        integer :: other_number = 0
        real(DP), dimension(Ndim) :: xCol
        real(DP), dimension(Ndim) :: mCol
        integer :: same_iCell
        integer :: mix_iCell
        logical :: add
    end type particle_index
    
    type, public :: particle_energy
        real(DP) :: same
        real(DP) :: mix
    end type particle_energy

    ! For main program

    type, public :: argument_random
        character(len=1) :: choice
        integer :: size
        integer, dimension(:), allocatable :: seed
    end type argument_random

    type, public :: argument_initial
        character(len=1) :: choice
        character(len=4096), dimension(3) :: files
        integer, dimension(3) :: length
    end type argument_initial
    
    type, public :: monteCarlo_arguments
        type(argument_random) :: random
        type(argument_initial) :: initial
    end type monteCarlo_arguments

end module module_types
