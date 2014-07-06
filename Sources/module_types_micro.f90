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
    
    ! Observables

    type, public :: Discrete_Observables
        integer :: num_hits = 0
        integer :: num_rejections = 0
        real(DP) :: rejection_rate = 0._DP
        real(DP) :: sum_rejection = 0._DP
    contains
        procedure :: update_rejection => Discrete_Observables_update_rejection
    end type Discrete_Observables
    
    type, extends(Discrete_Observables), public :: Adapting_Discrete_Observables
        real(DP) :: rejection_adapt = 0._DP
        real(DP) :: rejection_average = 0._DP
    contains
        procedure :: accumulate_rejection => Adapting_Discrete_Observables_accumulate_rejection
        procedure :: average_rejection => Adapting_Discrete_Observables_average_rejection
    end type Adapting_Discrete_Observables
    
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
    
contains

    subroutine Discrete_Observables_update_rejection(this)
    
        class(Discrete_Observables), intent(inout) :: this    

        this%rejection_rate = real(this%num_rejections, DP) / real(this%num_hits, DP)
        this%num_rejections = 0
        this%num_hits = 0
        
    end subroutine Discrete_Observables_update_rejection
    
    subroutine Adapting_Discrete_Observables_accumulate_rejection(this)
    
        class(Adapting_Discrete_Observables), intent(inout) :: this   
    
        this%rejection_adapt = this%rejection_adapt + this%rejection_rate
        
    end subroutine Adapting_Discrete_Observables_accumulate_rejection
    
    subroutine Adapting_Discrete_Observables_average_rejection(this, period_adaptation)
    
        class(Adapting_Discrete_Observables), intent(inout) :: this
        integer, intent(in) :: period_adaptation
    
        this%rejection_average = this%rejection_adapt / real(period_adaptation - 1, DP)
        this%rejection_adapt = 0._DP
        
    end subroutine Adapting_Discrete_Observables_average_rejection

end module module_types_micro
