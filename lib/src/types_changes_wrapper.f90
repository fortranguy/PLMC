module types_changes_wrapper

use data_constants, only: num_components
use class_moved_positions, only: Abstract_Moved_Positions
use class_rotated_orientations, only: Abstract_Rotated_Orientations
use class_change_tuner, only: Abstract_Change_Tuner
use class_particles_exchange, only: Abstract_Particles_Exchange

implicit none

private

    type, public :: Changes_Wrapper
        class(Abstract_Moved_Positions), allocatable :: moved_positions
        class(Abstract_Rotated_Orientations), allocatable :: rotated_orientations
        class(Abstract_Change_Tuner), allocatable :: move_tuner, rotation_tuner
        class(Abstract_Particles_Exchange), allocatable :: particles_exchange
    end type Changes_Wrapper

    type, public :: Mixture_Changes_Wrapper
        type(Changes_Wrapper) :: changes(num_components)
    end type Mixture_Changes_Wrapper

end module types_changes_wrapper
