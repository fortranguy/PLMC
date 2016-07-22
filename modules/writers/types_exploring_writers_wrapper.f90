module types_exploring_writers_wrapper

use classes_real_writer, only: Abstract_Real_Writer
use classes_line_writer, only: Abstract_Line_Writer
use classes_triangle_writer, only: Abstract_Triangle_Writer

implicit none

private

    type, public :: Exploring_Writers_Wrapper
        class(Abstract_Real_Writer), allocatable :: beta_pressure_excess
        class(Abstract_Line_Writer), allocatable :: inv_pow_activities
        class(Abstract_Line_Writer), allocatable :: field, walls
        class(Abstract_Triangle_Writer), allocatable :: short_energies, dipolar_energies
        class(Abstract_Real_Writer), allocatable :: dipolar_mixture_energy
        class(Abstract_Line_Writer), allocatable :: insertion_successes
    end type Exploring_Writers_Wrapper

end module types_exploring_writers_wrapper
