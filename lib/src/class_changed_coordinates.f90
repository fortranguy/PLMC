module class_changed_coordinates

implicit none

private

    type, abstract, public :: Abstract_Changed_Coordinates
    contains
        procedure(Abstract_Changed_Coordinates_increase_delta), deferred :: increase_delta
        procedure(Abstract_Changed_Coordinates_decrease_delta), deferred :: decrease_delta
    end type

    abstract interface

        subroutine Abstract_Changed_Coordinates_increase_delta(this)
        import :: Abstract_Changed_Coordinates
            class(Abstract_Changed_Coordinates), intent(inout) :: this
        end subroutine Abstract_Changed_Coordinates_increase_delta

        subroutine Abstract_Changed_Coordinates_decrease_delta(this)
        import :: Abstract_Changed_Coordinates
            class(Abstract_Changed_Coordinates), intent(inout) :: this
        end subroutine Abstract_Changed_Coordinates_decrease_delta

    end interface

end module class_changed_coordinates
