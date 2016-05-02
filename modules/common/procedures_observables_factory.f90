module procedures_observables_factory

use types_component_wrapper, only: Component_Wrapper
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use procedures_generating_observables_factory, only: generating_observables_create => create, &
    generating_observables_destroy => destroy
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper

implicit none

private
public :: observables_create_generating, observables_destroy_generating

contains

    subroutine observables_create_generating(observables, components)
        type(Generating_Observables_Wrapper), intent(out) :: observables
        type(Component_Wrapper), intent(in) :: components(:)

        call generating_observables_create(observables, components)
    end subroutine observables_create_generating

    subroutine observables_destroy_generating(observables)
        type(Generating_Observables_Wrapper), intent(out) :: observables

        call generating_observables_destroy(observables)
    end subroutine observables_destroy_generating

    !subroutine observables_create_exploring(observables)
    !    type(Exploring_Observables_Wrapper), intent(out) :: observables
    !end subroutine observables_create_exploring

end module procedures_observables_factory
