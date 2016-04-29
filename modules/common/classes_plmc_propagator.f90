module classes_plmc_propagator

use classes_tower_sampler, only: Abstract_Tower_Sampler
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithm_Pointer, &
    Metropolis_Algorithms_Wrapper

implicit none

private

    type, abstract, public :: Abstract_PLMC_Propagator
    private
        type(Metropolis_Algorithm_Pointer), pointer :: metropolis_algorithms(:) => null()
        class(Abstract_Tower_Sampler), allocatable :: selector
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
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

end module classes_plmc_propagator
