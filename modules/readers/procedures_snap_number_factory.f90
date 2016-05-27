module procedures_snap_number_factory

use classes_snap_number, only: Abstract_Snap_Number, Concrete_Snap_Number, Null_Snap_Number

implicit none

private
public :: create, destroy

contains

    subroutine create(snap_number, num_snaps, num_offset)
        class(Abstract_Snap_Number), allocatable, intent(out) :: snap_number
        integer, intent(in) :: num_snaps, num_offset

        if (num_snaps > 0) then
            allocate(Concrete_Snap_Number :: snap_number)
        else
            allocate(Null_Snap_Number :: snap_number)
        end if
        call snap_number%set(num_snaps, num_offset)
    end subroutine create

    subroutine destroy(snap_number)
        class(Abstract_Snap_Number), allocatable, intent(inout) :: snap_number

        if (allocated(snap_number)) deallocate(snap_number)
    end subroutine destroy

end module procedures_snap_number_factory
