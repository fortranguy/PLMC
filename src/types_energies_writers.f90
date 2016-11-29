module types_energies_writers

use classes_real_writer, only: Abstract_Real_Writer
use classes_line_writer, only: Abstract_Line_Writer
use classes_triangle_writer, only: Abstract_Triangle_Writer

implicit none

private

    type, public :: Concrete_Energies_Writers
        class(Abstract_Line_Writer), allocatable :: walls_energies, field_energies
        class(Abstract_Triangle_Writer), allocatable :: short_energies, dipolar_energies
        class(Abstract_Real_Writer), allocatable :: dipolar_shared_energy
    end type Concrete_Energies_Writers

end module types_energies_writers
