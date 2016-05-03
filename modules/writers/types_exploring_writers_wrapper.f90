module types_exploring_writers_wrapper

use classes_real_writer, only: Abstract_Real_Writer
use classes_line_writer, only: Abstract_Line_Writer

implicit none

private

    type, public :: Exploring_Writers_Wrapper
        class(Abstract_Real_Writer), allocatable :: pressure
        class(Abstract_Line_Writer), allocatable :: inv_pow_activities
    end type Exploring_Writers_Wrapper

end module types_exploring_writers_wrapper
