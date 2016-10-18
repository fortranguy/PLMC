module classes_generating_algorithm

use types_generating_observables_wrapper, only: Generating_Observables_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Generating_Algorithm
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_reset_selectors), deferred :: reset_selectors
        procedure(Abstract_get_num_choices), deferred :: get_num_choices
        procedure(Abstract_try), deferred :: try
    end type Abstract_Generating_Algorithm

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Generating_Algorithm
            class(Abstract_Generating_Algorithm), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_reset_selectors(this)
        import :: Abstract_Generating_Algorithm
            class(Abstract_Generating_Algorithm), intent(inout) :: this
        end subroutine Abstract_reset_selectors

        pure integer function Abstract_get_num_choices(this) result(num_choices)
        import :: Abstract_Generating_Algorithm
            class(Abstract_Generating_Algorithm), intent(in) :: this
        end function Abstract_get_num_choices

        subroutine Abstract_try(this, observables)
        import :: Generating_Observables_Wrapper, Abstract_Generating_Algorithm
            class(Abstract_Generating_Algorithm), intent(in) :: this
            type(Generating_Observables_Wrapper), intent(inout) :: observables
        end subroutine Abstract_try

    end interface

    type, extends(Abstract_Generating_Algorithm), public :: Null_Generating_Algorithm
    contains
        procedure :: destroy => Null_destroy
        procedure :: reset_selectors => Null_reset_selectors
        procedure :: get_num_choices => Null_get_num_choices
        procedure :: try => Null_try
    end type Null_Generating_Algorithm

    type, public :: Generating_Algorithm_Wrapper
        class(Abstract_Generating_Algorithm), allocatable :: algorithm
    end type Generating_Algorithm_Wrapper

contains

!implementation Null_Generating_Algorithm

    subroutine Null_destroy(this)
        class(Null_Generating_Algorithm), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_reset_selectors(this)
        class(Null_Generating_Algorithm), intent(inout) :: this
    end subroutine Null_reset_selectors

    pure integer function Null_get_num_choices(this) result(num_choices)
        class(Null_Generating_Algorithm), intent(in) :: this
        num_choices = 0
    end function Null_get_num_choices

    subroutine Null_try(this, observables)
        class(Null_Generating_Algorithm), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables
    end subroutine Null_try

!end implementation Null_Generating_Algorithm

end module classes_generating_algorithm
