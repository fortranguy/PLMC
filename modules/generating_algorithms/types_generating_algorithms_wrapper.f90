module types_generating_algorithms_wrapper

use classes_generating_algorithm, only: Abstract_Generating_Algorithm
use classes_box_volume_change, only: Abstract_Box_Volume_Change
use classes_box_particle_move, only: Abstract_Box_Particle_Move
use classes_box_particle_exchange, only: Box_Particle_Exchange
use classes_box_particles_swap, only: Abstract_Box_Particles_Swap

implicit none

private

    type, public :: Generating_Algorithm_Pointer
        class(Abstract_Generating_Algorithm), pointer :: algorithm => null()
    end type Generating_Algorithm_Pointer

    type, public :: Generating_Algorithms_Wrapper
        class(Abstract_Box_Volume_Change), allocatable :: box_volume_change
        class(Abstract_Generating_Algorithm), allocatable :: one_particle_translation, &
            one_particle_rotation
        class(Abstract_Box_Particles_Swap), allocatable :: two_particles_switch
        class(Abstract_Generating_Algorithm), allocatable :: one_particle_add, one_particle_remove
        class(Abstract_Box_Particles_Swap), allocatable :: two_particles_transmutation
    end type Generating_Algorithms_Wrapper

end module types_generating_algorithms_wrapper
