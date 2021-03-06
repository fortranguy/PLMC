module types_changes_wrapper

use classes_changed_box_size, only: Abstract_Changed_Box_Size
use classes_exchanged_boxes_size, only: Exchanged_Boxes_Size_Line
use classes_move_tuner, only: Abstract_Move_Tuner, Move_Tuner_Line
use classes_random_coordinates, only: Abstract_Random_Coordinates
use classes_coordinates_copier, only: Abstract_Coordinates_Copier
use types_changes_component_wrapper, only: Changes_Component_Wrapper

implicit none

private

    type, public :: Changes_Wrapper
        class(Abstract_Changed_Box_Size), allocatable :: changed_boxes_size(:)
        class(Abstract_Move_Tuner), allocatable :: boxes_size_change_tuner(:)
        type(Exchanged_Boxes_Size_Line), allocatable :: exchanged_boxes_size(:)
        type(Move_Tuner_Line), allocatable :: boxes_size_exchange_tuner(:)
        class(Abstract_Random_Coordinates), allocatable :: random_positions(:), random_orientation
        class(Abstract_Coordinates_Copier), allocatable :: position_copiers(:), orientation_copier
        type(Changes_Component_Wrapper), allocatable :: components(:, :)
    end type Changes_Wrapper

end module types_changes_wrapper
