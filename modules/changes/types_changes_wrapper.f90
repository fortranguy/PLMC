module types_changes_wrapper

use classes_changed_box_size, only: Changed_Box_Size_Line
use classes_move_tuner, only: Move_Tuner_Line
use classes_random_coordinates, only: Abstract_Random_Coordinates
use classes_coordinates_copier, only: Abstract_Coordinates_Copier
use types_changes_component_wrapper, only: Changes_Component_Wrapper

implicit none

private

    type, public :: Changes_Wrapper
        type(Changed_Box_Size_Line), allocatable :: changed_boxes_size(:)
        type(Move_Tuner_Line), allocatable :: boxes_size_change_tuner(:)
        class(Abstract_Random_Coordinates), allocatable :: random_positions(:), random_orientation
        class(Abstract_Coordinates_Copier), allocatable :: position_copiers(:), orientation_copier
        type(Changes_Component_Wrapper), allocatable :: components(:, :)
    end type Changes_Wrapper

end module types_changes_wrapper
