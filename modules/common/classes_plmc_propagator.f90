module classes_plmc_propagator

use classes_tower_sampler, only: Concrete_Tower_Sampler
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithm_Pointer, &
    Metropolis_Algorithms_Wrapper
use types_observables_wrapper, only: Observables_Wrapper

implicit none

private

    type, abstract, public :: Abstract_PLMC_Propagator
    private
        type(Metropolis_Algorithm_Pointer), allocatable :: metropolis_algorithms(:)
        type(Concrete_Tower_Sampler) :: selector
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: set_selector => Abstract_set_selector
    end type Abstract_PLMC_Propagator

    type, extends(Abstract_PLMC_Propagator), public :: Concrete_PLMC_Propagator

    end type Concrete_PLMC_Propagator

contains

    subroutine Abstract_construct(this, metropolis_algorithms)
        class(Abstract_PLMC_Propagator), intent(out) :: this
        type(Metropolis_Algorithm_Pointer), target, intent(in) :: metropolis_algorithms(:)

        integer :: i_algorithm

        allocate(this%metropolis_algorithms(size(metropolis_algorithms))) !shorter with source?
        do i_algorithm = 1, size(this%metropolis_algorithms)
            this%metropolis_algorithms(i_algorithm)%algorithm => &
                metropolis_algorithms(i_algorithm)%algorithm
        end do
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_PLMC_Propagator), intent(inout) :: this

        integer :: i_algorithm

        call this%selector%destroy()
        if (allocated(this%metropolis_algorithms)) then
            do i_algorithm = size(this%metropolis_algorithms), 1, -1
                this%metropolis_algorithms(i_algorithm)%algorithm => null()
            end do
            deallocate(this%metropolis_algorithms)
        end if
    end subroutine Abstract_destroy

    subroutine Abstract_set_selector(this)
        class(Abstract_PLMC_Propagator), intent(inout) :: this

        integer :: nums_choices(size(this%metropolis_algorithms)), i_choice

        do i_choice = 1, size(nums_choices)
            nums_choices(i_choice) = this%metropolis_algorithms(i_choice)%algorithm%&
                get_num_choices()
        end do
        call this%selector%construct(nums_choices)
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
