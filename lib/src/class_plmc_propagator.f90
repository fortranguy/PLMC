module class_plmc_propagator

use procedures_checks, only: check_positive
use class_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use class_metropolis_algorithm, only: Abstract_Metropolis_Algorithm
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithm_Pointer
use types_observables_wrapper, only: Observables_Wrapper
implicit none

private

    type, abstract, public :: Abstract_PLMC_Propagator
    private
        type(Metropolis_Algorithm_Pointer), allocatable :: metropolis_algorithms(:)
        class(Abstract_Tower_Sampler), allocatable :: selector
        integer :: num_choices
    contains
        procedure :: construct => Abstract_PLMC_Propagator_construct
        procedure :: destroy => Abstract_PLMC_Propagator_destroy
        procedure :: try => Abstract_PLMC_Propagator_try
    end type Abstract_PLMC_Propagator

    type, extends(Abstract_PLMC_Propagator), public :: Concrete_PLMC_Propagator

    end type Concrete_PLMC_Propagator

contains

!implementation Abstract_PLMC_Propagator

    subroutine Abstract_PLMC_Propagator_construct(this, metropolis_algorithms)
        class(Abstract_PLMC_Propagator), intent(inout) :: this
        type(Metropolis_Algorithm_Pointer), intent(in) :: metropolis_algorithms(:)

        integer, allocatable :: nums_choices(:)
        integer :: i_choice

        call check_positive("Abstract_PLMC_Propagator_construct", "size(metropolis_algorithms)", &
            size(metropolis_algorithms))
        allocate(this%metropolis_algorithms(size(metropolis_algorithms)))
        this%metropolis_algorithms = metropolis_algorithms

        allocate(nums_choices(size(this%metropolis_algorithms)))
        do i_choice = 1, size(nums_choices)
            nums_choices(i_choice) = this%metropolis_algorithms(i_choice)%algorithm%&
                get_num_choices()
        end do
        this%num_choices = sum(nums_choices)
        if (this%num_choices == 0) then
            allocate(Null_Tower_Sampler :: this%selector)
        else
            allocate(Concrete_Tower_Sampler :: this%selector)
        end if
        call this%selector%construct(nums_choices)
        deallocate(nums_choices)
    end subroutine Abstract_PLMC_Propagator_construct

    subroutine Abstract_PLMC_Propagator_destroy(this)
        class(Abstract_PLMC_Propagator), intent(inout) :: this

        if (allocated(this%selector)) then
            call this%selector%destroy()
            deallocate(this%selector)
        end if
        if (allocated(this%metropolis_algorithms)) deallocate(this%metropolis_algorithms)
    end subroutine Abstract_PLMC_Propagator_destroy

    subroutine Abstract_PLMC_Propagator_try(this, observables)
        class(Abstract_PLMC_Propagator), intent(in) :: this
        type(Observables_Wrapper), intent(inout) :: observables

        integer :: i_choice

        do i_choice = 1, this%num_choices
            call this%metropolis_algorithms(this%selector%get())%algorithm%try(observables)
        end do
    end subroutine Abstract_PLMC_Propagator_try

!end implementation Abstract_PLMC_Propagator

end module class_plmc_propagator
