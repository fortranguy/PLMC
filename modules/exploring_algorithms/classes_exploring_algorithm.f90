module classes_exploring_algorithm

use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Exploring_Algorithm
    contains
        procedure(Abstract_try), deferred :: try
    end type Abstract_Exploring_Algorithm

    abstract interface

        subroutine Abstract_try(this, observables)
        import :: Exploring_Observables_Wrapper, Abstract_Exploring_Algorithm
            class(Abstract_Exploring_Algorithm), intent(in) :: this
            type(Exploring_Observables_Wrapper), intent(inout) :: observables
        end subroutine Abstract_try

    end interface

end module classes_exploring_algorithm
