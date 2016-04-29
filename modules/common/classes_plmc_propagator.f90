module classes_plmc_propagator

use classes_tower_sampler, only: Abstract_Tower_Sampler
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithm_Pointer, &
    Metropolis_Algorithms_Wrapper
use types_observables_wrapper, only: Observables_Wrapper

implicit none

private

    type, abstract, public :: Abstract_PLMC_Propagator
    private
        type(Metropolis_Algorithm_Pointer), pointer :: metropolis_algorithms(:) => null()
        class(Abstract_Tower_Sampler), allocatable :: selector
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: set_selector => Abstract_set_selector
    end type Abstract_PLMC_Propagator

contains

    subroutine Abstract_construct(this, metropolis_algorithms)
        class(Abstract_PLMC_Propagator), intent(out) :: this
        type(Metropolis_Algorithm_Pointer), target, intent(in) :: metropolis_algorithms(:)

        this%metropolis_algorithms => metropolis_algorithms
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_PLMC_Propagator), intent(inout) :: this

        if (allocated(this%selector)) then
            call this%selector%destroy()
            deallocate(this%selector)
        end if
        this%metropolis_algorithms => null()
    end subroutine Abstract_destroy

    subroutine Abstract_set_selector(this, selector)
        class(Abstract_PLMC_Propagator), intent(inout) :: this
        class(Abstract_Tower_Sampler), intent(in) :: selector

        allocate(this%selector, source=selector)
    end subroutine Abstract_set_selector

    subroutine Abstract_try(this, observables)
        class(Abstract_PLMC_Propagator), intent(in) :: this
        type(Observables_Wrapper), intent(inout) :: observables

        integer :: i_choice, i_random

        do i_choice = 1, this%selector%get_num_choices()
            i_random = this%selector%get() !gfortran bug?
            call this%metropolis_algorithms(i_random)%algorithm%try(observables)
        end do
    end subroutine Abstract_try

end module classes_plmc_propagator
