module types_generating_writers_wrapper

use classes_box_size_writer, only: Abstract_Box_Size_Writer
use classes_real_writer, only: Abstract_Real_Writer
use classes_line_writer, only: Abstract_Line_Writer
use classes_triangle_writer, only: Abstract_Triangle_Writer
use types_component_writers_wrapper, only: Component_Writers_Wrapper

implicit none

private

    type, public :: Generating_Writers_Wrapper
        class(Abstract_Box_Size_Writer), allocatable :: box_size
        class(Abstract_Line_Writer), allocatable :: field, walls
        type(Component_Writers_Wrapper), allocatable :: components(:)
        class(Abstract_Triangle_Writer), allocatable :: switches, short_energies, dipolar_energies
        class(Abstract_Real_Writer), allocatable :: dipolar_mixture_energy
    end type Generating_Writers_Wrapper

end module types_generating_writers_wrapper
