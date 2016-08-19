module classes_generating_algorithm

use types_generating_observables_wrapper, only: Generating_Observables_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Generating_Algorithm
    contains
        procedure(Abstract_get_num_choices), deferred :: get_num_choices
        procedure(Abstract_try), deferred :: try
    end type Abstract_Generating_Algorithm

    abstract interface

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

end module classes_generating_algorithm
