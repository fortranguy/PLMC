module types_move_tuner_wrapper

use classes_move_tuner, only: Abstract_Move_Tuner

implicit none

private

    type, public :: Move_Tuner_Wrapper
        class(Abstract_Move_Tuner), allocatable :: tuner
    end type Move_Tuner_Wrapper

    type, public :: Move_Tuner_Line
        type(Move_Tuner_Wrapper), allocatable :: line(:)
    end type Move_Tuner_Line

end module types_move_tuner_wrapper
