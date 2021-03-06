module types_generating_writers_wrapper

use classes_real_writer, only: Abstract_Real_Writer
use classes_line_writer, only: Abstract_Line_Writer, Line_Writer_Wrapper
use classes_triangle_writer, only: Abstract_Triangle_Writer
use classes_rectangle_writer, only: Abstract_Rectangle_Writer, Rectangle_Writer_Wrapper
use classes_complete_coordinates_writer, only: Abstract_Complete_Coordinates_Writer
use classes_changes_success_writer, only: Changes_Success_Writer_Wrapper
use types_energies_writers, only: Concrete_Energies_Writers

implicit none

private

    type, public :: Generating_Writers_Wrapper
        class(Abstract_Real_Writer), allocatable :: accessible_domains_size(:)
        class(Abstract_Line_Writer), allocatable :: volumes_change_success
        class(Abstract_Triangle_Writer), allocatable :: volumes_exchange_success
        type(Line_Writer_Wrapper), allocatable :: teleportations_successes(:, :)
        type(Rectangle_Writer_Wrapper), allocatable :: swaps_successes(:, :)
        class(Abstract_Line_Writer), allocatable :: nums_particles(:)
        class(Abstract_Complete_Coordinates_Writer), allocatable :: complete_coordinates
        type(Concrete_Energies_Writers), allocatable :: energies(:)
        type(Changes_Success_Writer_Wrapper), allocatable :: components_changes(:, :)
        class(Abstract_Triangle_Writer), allocatable :: switches_successes(:)
        class(Abstract_Rectangle_Writer), allocatable :: transmutations_successes(:)
    end type Generating_Writers_Wrapper

end module types_generating_writers_wrapper
