module types_metropolis_algorithms_wrapper

use class_metropolis_algorithm, only: Abstract_Metropolis_Algorithm
use class_one_particle_change, only: Abstract_One_Particle_Change

implicit none

private

    type, public :: Metropolis_Algorithm_Pointer
        class(Abstract_Metropolis_Algorithm), pointer :: algorithm => null()
    end type Metropolis_Algorithm_Pointer

    type, public :: Metropolis_Algorithms_Wrapper
        class(Abstract_One_Particle_Change), allocatable :: one_particle_move, one_particle_rotation
    end type Metropolis_Algorithms_Wrapper

end module types_metropolis_algorithms_wrapper
