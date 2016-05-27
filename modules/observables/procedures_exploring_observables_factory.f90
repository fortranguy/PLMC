module procedures_exploring_observables_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_component_wrapper, only: Component_Wrapper
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_all
    module procedure :: create_inv_pow_activities
end interface create

interface destroy
    module procedure :: destroy_inv_pow_activities
    module procedure :: destroy_all
end interface destroy

contains

    subroutine create_all(observables, components)
        type(Exploring_Observables_Wrapper), intent(out) ::observables
        type(Component_Wrapper), intent(in) :: components(:)

        call create(observables%inv_pow_activities, size(components))
    end subroutine create_all

    subroutine destroy_all(observables)
        type(Exploring_Observables_Wrapper), intent(inout) ::observables

        call destroy(observables%inv_pow_activities)
    end subroutine destroy_all

    subroutine create_inv_pow_activities(inv_pow_activities, num_components)
        real(DP), allocatable, intent(out) :: inv_pow_activities(:)
        integer, intent(in) :: num_components

        allocate(inv_pow_activities(num_components))
    end subroutine create_inv_pow_activities

    subroutine destroy_inv_pow_activities(inv_pow_activities)
        real(DP), allocatable, intent(inout) :: inv_pow_activities(:)

        if (allocated(inv_pow_activities)) deallocate(inv_pow_activities)
    end subroutine destroy_inv_pow_activities

end module procedures_exploring_observables_factory
