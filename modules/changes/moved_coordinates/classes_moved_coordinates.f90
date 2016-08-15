module classes_moved_coordinates

implicit none

private

    type, abstract, public :: Abstract_Moved_Coordinates
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_increase_delta), deferred :: increase_delta
        procedure(Abstract_decrease_delta), deferred :: decrease_delta
    end type Abstract_Moved_Coordinates

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Moved_Coordinates
            class(Abstract_Moved_Coordinates), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_increase_delta(this)
        import :: Abstract_Moved_Coordinates
            class(Abstract_Moved_Coordinates), intent(inout) :: this
        end subroutine Abstract_increase_delta

        subroutine Abstract_decrease_delta(this)
        import :: Abstract_Moved_Coordinates
            class(Abstract_Moved_Coordinates), intent(inout) :: this
        end subroutine Abstract_decrease_delta

    end interface

end module classes_moved_coordinates
