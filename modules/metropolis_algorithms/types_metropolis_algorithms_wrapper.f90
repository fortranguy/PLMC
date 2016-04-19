module types_metropolis_algorithms_wrapper

use classes_metropolis_algorithm, only: Abstract_Metropolis_Algorithm
use classes_one_particle_change, only: Abstract_One_Particle_Change
use classes_two_particles_switch, only: Abstract_Two_Particles_Switch

implicit none

private

    type, public :: Metropolis_Algorithm_Pointer
        class(Abstract_Metropolis_Algorithm), pointer :: algorithm => null()
    end type Metropolis_Algorithm_Pointer

    type, public :: Metropolis_Algorithms_Wrapper
        class(Abstract_One_Particle_Change), allocatable :: one_particle_move, one_particle_rotation
        class(Abstract_Two_Particles_Switch), allocatable :: two_particles_switch
    end type Metropolis_Algorithms_Wrapper

end module types_metropolis_algorithms_wrapper
