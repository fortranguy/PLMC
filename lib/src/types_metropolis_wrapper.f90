module types_metropolis_wrapper

use class_one_particle_change, only: Abstract_One_Particle_Change

implicit none

private

    integer, public, parameter :: num_algorithms = 2

    type, public :: Metropolis_Wrapper
        class(Abstract_One_Particle_Change), allocatable :: one_particle_move, one_particle_rotation
    end type Metropolis_Wrapper

end module types_metropolis_wrapper
