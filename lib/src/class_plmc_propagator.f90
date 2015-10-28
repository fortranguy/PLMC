module class_plmc_propagator

use procedures_errors, only: error_exit
use class_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use class_metropolis_algorithm, only: Abstract_Metropolis_Algorithm
use types_metropolis_wrapper, only: Metropolis_Algorithm_Pointer
use types_observables_wrapper, only: Observables_Wrapper
implicit none

private

    type, abstract, public :: Abstract_PLMC_Propagator
    private
        type(Metropolis_Algorithm_Pointer), allocatable :: metropolis(:)
        class(Abstract_Tower_Sampler), allocatable :: selector
        integer :: num_choices
    contains
        procedure :: add => Abstract_PLMC_Propagator_add
        procedure :: construct => Abstract_PLMC_Propagator_construct
        procedure :: destroy => Abstract_PLMC_Propagator_destroy
        procedure :: try => Abstract_PLMC_Propagator_try
    end type Abstract_PLMC_Propagator

    type, extends(Abstract_PLMC_Propagator), public :: Concrete_PLMC_Propagator

    end type Concrete_PLMC_Propagator

contains

!implementation Abstract_PLMC_Propagator

    subroutine Abstract_PLMC_Propagator_add(this, algorithm)
        class(Abstract_PLMC_Propagator), intent(inout) :: this
        class(Abstract_Metropolis_Algorithm), target, intent(in) :: algorithm

        type(Metropolis_Algorithm_Pointer) :: new

        if (.not.allocated(this%metropolis)) then
            allocate(this%metropolis(1))
            this%metropolis(1)%algorithm => algorithm
        else
            new%algorithm => algorithm
            this%metropolis = [this%metropolis, new]
        end if
    end subroutine Abstract_PLMC_Propagator_add

    subroutine Abstract_PLMC_Propagator_construct(this)
        class(Abstract_PLMC_Propagator), intent(inout) :: this

        integer :: nums_choices(size(this%metropolis)), i_choice

        do i_choice = 1, size(this%metropolis)
            nums_choices(i_choice) = this%metropolis(i_choice)%algorithm%get_num_choices()
        end do
        this%num_choices = sum(nums_choices)
        if (this%num_choices == 0) then
            allocate(Null_Tower_Sampler :: this%selector)
        else
            allocate(Concrete_Tower_Sampler :: this%selector)
        end if
        call this%selector%construct(nums_choices)
    end subroutine Abstract_PLMC_Propagator_construct

    subroutine Abstract_PLMC_Propagator_destroy(this)
        class(Abstract_PLMC_Propagator), intent(inout) :: this

        if (allocated(this%selector)) then
            call this%selector%destroy()
            deallocate(this%selector)
        end if
    end subroutine Abstract_PLMC_Propagator_destroy

    subroutine Abstract_PLMC_Propagator_try(this, observables)
        class(Abstract_PLMC_Propagator), intent(in) :: this
        type(Observables_Wrapper), intent(inout) :: observables

        integer :: i_choice

        do i_choice = 1, this%num_choices
            call this%metropolis(this%selector%get())%algorithm%try(observables)
        end do
    end subroutine Abstract_PLMC_Propagator_try

!end implementation Abstract_PLMC_Propagator

end module class_plmc_propagator
