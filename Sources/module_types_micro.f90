!> \brief Definition of derived types (micro)

module module_types_micro

use data_precisions, only: DP
use data_box, only: num_dimensions

implicit none

    ! System Box
    
    type, public :: Box_Parameters
        real(DP), dimension(num_dimensions) :: size
        integer, dimension(num_dimensions) :: wave
        real(DP) :: temperature
        integer :: num_particles
    end type Box_Parameters
    
    ! Algorithms
    
    type, public :: Particle_Index
        integer :: number
        integer :: other_number = 0
        real(DP), dimension(num_dimensions) :: position
        real(DP), dimension(num_dimensions) :: orientation
        integer :: same_i_cell
        integer :: between_i_cell
        logical :: add
    end type Particle_Index
    
    type, public :: Particle_Energy
        real(DP) :: same
        real(DP) :: mix
    end type Particle_Energy
    
    ! Neighbour cells

    type, public :: Node
        integer :: number
        type(Node), pointer :: next => null()
    end type Node

    type, public :: Linked_List
        type(Node), pointer :: particle => null()
    end type Linked_List

    ! Arguments

    type, public :: Argument_Random
        character(len=1) :: choice
        integer :: size
        integer, dimension(:), allocatable :: seed
    end type Argument_Random

    type, public :: Argument_Initial
        character(len=1) :: choice
        character(len=4096), dimension(3) :: files
        integer, dimension(3) :: length
    end type Argument_Initial
    
    type, public :: Monte_Carlo_Arguments
        type(Argument_Random) :: random
        type(Argument_Initial) :: initial
    end type Monte_Carlo_Arguments

end module module_types_micro
