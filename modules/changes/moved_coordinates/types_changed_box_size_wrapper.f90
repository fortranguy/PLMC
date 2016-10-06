module types_changed_box_size_wrapper

use classes_changed_box_size, only: Abstract_Changed_Box_Size

implicit none

private

    type, public :: Changed_Box_Size_Wrapper
        class(Abstract_Changed_Box_Size), allocatable :: changed
    end type Changed_Box_Size_Wrapper

    type, public :: Changed_Box_Size_Line
        class(Changed_Box_Size_Wrapper), allocatable :: line(:)
    end type Changed_Box_Size_Line

end module types_changed_box_size_wrapper
