module types_generating_algorithms_wrapper

use classes_generating_algorithm, only: Abstract_Generating_Algorithm
use classes_box_volume_change, only: Abstract_Box_Volume_Change
use classes_one_particle_move, only: Abstract_One_Particle_Move
use classes_two_particles_switch, only: Abstract_Two_Particles_Switch
use classes_one_particle_exchange, only: Abstract_One_Particle_Exchange
use classes_two_particles_transmutation, only: Abstract_Two_Particles_Transmutation

implicit none

private

    type, public :: Generating_Algorithm_Pointer
        class(Abstract_Generating_Algorithm), pointer :: algorithm => null()
    end type Generating_Algorithm_Pointer

    type, public :: Generating_Algorithms_Wrapper
        class(Abstract_Box_Volume_Change), allocatable :: box_volume_change
        class(Abstract_One_Particle_Move), allocatable :: one_particle_translation, &
            one_particle_rotation
        class(Abstract_Two_Particles_Switch), allocatable :: two_particles_switch
        class(Abstract_One_Particle_Exchange), allocatable :: one_particle_add, one_particle_remove
        class(Abstract_Two_Particles_Transmutation), allocatable :: two_particles_transmutation
    end type Generating_Algorithms_Wrapper

end module types_generating_algorithms_wrapper
