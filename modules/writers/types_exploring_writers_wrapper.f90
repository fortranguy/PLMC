module types_exploring_writers_wrapper

use classes_real_writer, only: Abstract_Real_Writer
use classes_line_writer, only: Abstract_Line_Writer
use types_energies_writers, only: Concrete_Energies_Writers

implicit none

private

    type, public :: Exploring_Writers_Wrapper
        class(Abstract_Real_Writer), allocatable :: maximum_boxes_compression_delta(:)
        class(Abstract_Real_Writer), allocatable :: beta_pressures_excess(:)
        type(Concrete_Energies_Writers), allocatable :: energies(:)
        class(Abstract_Line_Writer), allocatable :: inv_pow_activities(:)
        class(Abstract_Line_Writer), allocatable :: insertion_successes(:)
    end type Exploring_Writers_Wrapper

end module types_exploring_writers_wrapper
