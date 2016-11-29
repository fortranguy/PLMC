module procedures_observables_factory

use types_component_wrapper, only: Component_Wrapper
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use procedures_generating_observables_factory, only: generating_create => create, &
    generating_destroy => destroy
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use procedures_exploring_observables_factory, only: exploring_create => create, &
    exploring_destroy => destroy

implicit none

private
public :: observables_create_generating, observables_create_exploring, &
    observables_destroy_generating, observables_destroy_exploring

contains

    subroutine observables_create_generating(observables,  components)
        type(Generating_Observables_Wrapper), intent(out) :: observables
        type(Component_Wrapper), intent(in) :: components(:, :)

        call generating_create(observables, size(components, 2), size(components, 1))
    end subroutine observables_create_generating

    subroutine observables_destroy_generating(observables)
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        call generating_destroy(observables)
    end subroutine observables_destroy_generating

    subroutine observables_create_exploring(observables, components)
        type(Exploring_Observables_Wrapper), intent(out) :: observables
        type(Component_Wrapper), intent(in) :: components(:, :)

        call exploring_create(observables, size(components, 2), size(components, 1))
    end subroutine observables_create_exploring

    subroutine observables_destroy_exploring(observables)
        type(Exploring_Observables_Wrapper), intent(inout) :: observables

        call exploring_destroy(observables)
    end subroutine observables_destroy_exploring

end module procedures_observables_factory
