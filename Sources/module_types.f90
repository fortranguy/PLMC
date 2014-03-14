!> \brief Definition of derives types

module module_types

implicit none

    ! For neighbour cells

    type, public :: Node
        integer :: iCol
        type(Node), pointer :: next => null()
    end type Node

    type, public :: LinkedList
        type(Node), pointer :: particle => null()
    end type LinkedList

    ! For main program

    type, public :: argument_seed
        character(len=1) :: choice
        integer :: size
        integer, dimension(:), allocatable :: seed
    end type argument_seed

    type, public :: argument_initial
        character(len=1) :: choice
        character(len=4096), dimension(3) :: files
        integer, dimension(3) :: length
    end type argument_initial
    
    type, public :: monteCarlo_arguments
        type(argument_seed) :: seed
        type(argument_initial) :: initial
    end type monteCarlo_arguments

end module module_types
