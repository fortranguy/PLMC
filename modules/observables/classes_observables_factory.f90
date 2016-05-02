module classes_observables_factory

use types_observables_wrapper, only: Generating_Observables_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Observables_Factory
    contains
        procedure(Abstract_create), deferred, nopass :: create
        procedure(Abstract_destroy), deferred, nopass :: destroy
    end type Abstract_Observables_Factory

    abstract interface

        subroutine Abstract_create(observables)
        import :: Generating_Observables_Wrapper
            type(Generating_Observables_Wrapper), intent(out) :: observables
        end subroutine Abstract_create

        subroutine Abstract_destroy(observables)
        import :: Generating_Observables_Wrapper
            type(Generating_Observables_Wrapper), intent(inout) :: observables
        end subroutine Abstract_destroy

    end interface

end module classes_observables_factory
