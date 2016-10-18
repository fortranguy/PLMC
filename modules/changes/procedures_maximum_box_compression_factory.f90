module procedures_maximum_box_compression_factory

use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_environment_inquirers, only: periodicity_is_xyz, periodicity_is_xy
use classes_maximum_box_compression, only: Abstract_Maximum_Box_Compression, &
    XYZ_Maximum_Box_Compression, XY_Maximum_Box_Compression, Null_Maximum_Box_Compression

implicit none

private
public :: create, destroy

contains

    subroutine create(maximum_box_compression, periodic_boxes, measure)
        class(Abstract_Maximum_Box_Compression), allocatable, intent(out) :: maximum_box_compression
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        logical, intent(in) :: measure

        if (measure) then
            if (all(periodicity_is_xyz(periodic_boxes))) then
                allocate(XYZ_Maximum_Box_Compression :: maximum_box_compression)
            else if (all(periodicity_is_xy(periodic_boxes))) then
                allocate(XY_Maximum_Box_Compression :: maximum_box_compression)
            else
                call error_exit("procedures_maximum_box_compression_factory: create: "//&
                        "box periodicity is unknown.")
            end if
        else
            allocate(Null_Maximum_Box_Compression :: maximum_box_compression)
        end if
    end subroutine create

    subroutine destroy(maximum_box_compression)
        class(Abstract_Maximum_Box_Compression), allocatable, intent(inout) :: &
            maximum_box_compression

        if (allocated(maximum_box_compression)) deallocate(maximum_box_compression)
    end subroutine destroy

end module procedures_maximum_box_compression_factory
