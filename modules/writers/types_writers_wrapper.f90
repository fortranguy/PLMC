module types_writers_wrapper

use class_component_coordinates_writer, only: Abstract_Coordinates_Writer
use class_line_writer, only: Abstract_Line_Writer
use class_triangle_writer, only: Abstract_Triangle_Writer
use class_energy_writer, only: Abstract_Energy_Writer
use class_changes_writer, only: Abstract_Changes_Success_Writer

implicit none

private

    type, public :: Component_Writers_Wrapper
        class(Abstract_Coordinates_Writer), allocatable :: coordinates
        class(Abstract_Changes_Success_Writer), allocatable :: changes
    end type Component_Writers_Wrapper

    type, public :: Writers_Wrapper
        type(Component_Writers_Wrapper), allocatable :: components(:)
        class(Abstract_Line_Writer), allocatable :: walls
        class(Abstract_Triangle_Writer), allocatable :: switches, short_energies, long_energies
        class(Abstract_Energy_Writer), allocatable :: long_mixture_energy
    end type Writers_Wrapper

end module types_writers_wrapper
