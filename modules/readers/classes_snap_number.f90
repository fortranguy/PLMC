module classes_snap_number

implicit none

private

    type, abstract, public :: Abstract_Snap_Number
    private
        integer :: num_snaps
        integer :: num_offset
    contains
        procedure :: set => Abstract_set
        procedure :: get => Abstract_get
    end type Abstract_Snap_Number

    type, extends(Abstract_Snap_Number), public :: Concrete_Snap_Number

    end type Concrete_Snap_Number

    type, extends(Abstract_Snap_Number), public :: Null_Snap_Number
    contains
        procedure :: set => Null_set
        procedure :: get => Null_get
    end type Null_Snap_Number

contains

!implementation Abstract_Snap_Number

    subroutine Abstract_set(this, num_snaps, num_offset)
        class(Abstract_Snap_Number), intent(inout) :: this
        integer, intent(in) :: num_snaps, num_offset

        this%num_snaps = num_snaps
        this%num_offset = num_offset
    end subroutine Abstract_set

    pure integer function Abstract_get(this, i_component, i_snap) result(snap_number)
        class(Abstract_Snap_Number), intent(in) :: this
        integer, intent(in) :: i_component, i_snap

        snap_number = (i_component - 1)*this%num_snaps + i_snap + this%num_offset
    end function Abstract_get

!end implementation Abstract_Snap_Number

!implementation Null_Snap_Number

    subroutine Null_set(this, num_snaps, num_offset)
        class(Null_Snap_Number), intent(inout) :: this
        integer, intent(in) :: num_snaps, num_offset
    end subroutine Null_set

    pure integer function Null_get(this, i_component, i_snap) result(snap_number)
        class(Null_Snap_Number), intent(in) :: this
        integer, intent(in) :: i_component, i_snap
        snap_number = 0
    end function Null_get

!end implementation Null_Snap_Number

end module classes_snap_number
