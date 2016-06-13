module types_metropolis_algorithms_wrapper

use classes_metropolis_algorithm, only: Abstract_Metropolis_Algorithm
use classes_one_particle_move, only: Abstract_One_Particle_Move
use classes_two_particles_switch, only: Abstract_Two_Particles_Switch
use classes_one_particle_exchange, only: Abstract_One_Particle_Exchange
use classes_two_particles_transmutation, only: Abstract_Two_Particles_Transmutation

implicit none

private

    type, public :: Metropolis_Algorithm_Pointer
        class(Abstract_Metropolis_Algorithm), pointer :: algorithm => null()
    end type Metropolis_Algorithm_Pointer

    type, public :: Metropolis_Algorithms_Wrapper
        class(Abstract_One_Particle_Move), allocatable :: one_particle_translation, &
            one_particle_rotation
        class(Abstract_Two_Particles_Switch), allocatable :: two_particles_switch
        class(Abstract_One_Particle_Exchange), allocatable :: one_particle_add, one_particle_remove
        class(Abstract_Two_Particles_Transmutation), allocatable :: two_particles_transmutation
    end type Metropolis_Algorithms_Wrapper

end module types_metropolis_algorithms_wrapper
