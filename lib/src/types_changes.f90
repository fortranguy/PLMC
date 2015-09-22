module types_changes

use class_moved_positions, only: Abstract_Moved_Positions
use class_rotated_orientations, only: Abstract_Rotated_Orientations
use class_particles_exchange, only: Abstract_Particles_Exchange

implicit none

private

    type, public :: Changes_Wrapper
        class(Abstract_Moved_Positions), allocatable :: moved_positions
        class(Abstract_Rotated_Orientations), allocatable :: rotated_orientations
        class(Abstract_Particles_Exchange), allocatable :: particles_exchange
    end type Changes_Wrapper

end module types_changes