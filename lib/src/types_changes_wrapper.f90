module types_changes_wrapper

use data_constants, only: num_components
use class_changed_coordinates, only: Abstract_Changed_Coordinates
use class_change_tuner, only: Abstract_Change_Tuner
use class_component_exchange, only: Abstract_Component_Exchange

implicit none

private

    type, public :: Changes_Wrapper
        class(Abstract_Changed_Coordinates), allocatable :: moved_positions, rotated_orientations
        class(Abstract_Change_Tuner), allocatable :: move_tuner, rotation_tuner
        class(Abstract_Component_Exchange), allocatable :: component_exchange
    end type Changes_Wrapper

end module types_changes_wrapper
