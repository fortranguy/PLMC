module class_metropolis_algorithm

use types_observables_wrapper, only: Observables_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Metropolis_Algorithm
    contains
        procedure(Abstract_get_num_choices), deferred :: get_num_choices
        procedure(Abstract_try), deferred :: try
    end type Abstract_Metropolis_Algorithm

    abstract interface

        pure integer function Abstract_get_num_choices(this) result(num_choices)
        import :: Abstract_Metropolis_Algorithm
            class(Abstract_Metropolis_Algorithm), intent(in) :: this
        end function Abstract_get_num_choices

        subroutine Abstract_try(this, observables)
        import :: Observables_Wrapper, Abstract_Metropolis_Algorithm
            class(Abstract_Metropolis_Algorithm), intent(in) :: this
            type(Observables_Wrapper), intent(inout) :: observables
        end subroutine Abstract_try

    end interface

end module class_metropolis_algorithm
