module classes_tunable_move

implicit none

private

    type, abstract, public :: Abstract_Tunable_Move
    contains
        procedure(Abstract_increase_delta), deferred :: increase_delta
        procedure(Abstract_decrease_delta), deferred :: decrease_delta
    end type Abstract_Tunable_Move

    abstract interface

        subroutine Abstract_increase_delta(this)
        import :: Abstract_Tunable_Move
            class(Abstract_Tunable_Move), intent(inout) :: this
        end subroutine Abstract_increase_delta

        subroutine Abstract_decrease_delta(this)
        import :: Abstract_Tunable_Move
            class(Abstract_Tunable_Move), intent(inout) :: this
        end subroutine Abstract_decrease_delta

    end interface

end module classes_tunable_move
