module types_changes_component_wrapper

use classes_moved_component_coordinates, only: Abstract_Moved_Component_Coordinates
use classes_move_tuner, only: Abstract_Move_Tuner

implicit none

private

    type, public :: Changes_Component_Wrapper
        class(Abstract_Moved_Component_Coordinates), allocatable :: translated_positions, &
            rotated_orientations
        class(Abstract_Move_Tuner), allocatable :: translation_tuner, rotation_tuner
    end type Changes_Component_Wrapper

end module types_changes_component_wrapper
