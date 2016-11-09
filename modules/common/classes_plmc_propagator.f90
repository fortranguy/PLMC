module classes_plmc_propagator

use procedures_checks, only: check_positive
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only: tower_sampler_destroy => destroy
use classes_generating_algorithm, only: Generating_Algorithm_Wrapper
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper

implicit none

private

    type, abstract, public :: Abstract_PLMC_Propagator
    private
        type(Generating_Algorithm_Wrapper), pointer :: generating_algorithms(:) => null()
        class(Abstract_Tower_Sampler), allocatable :: selector
        integer :: tuning_period = 0
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset => Abstract_reset
        procedure :: tune => Abstract_tune
        procedure :: try => Abstract_try
    end type Abstract_PLMC_Propagator

    type, extends(Abstract_PLMC_Propagator), public :: Concrete_PLMC_Propagator

    end type Concrete_PLMC_Propagator

    type, extends(Abstract_PLMC_Propagator), public :: Null_PLMC_Propagator
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: reset => Null_reset
        procedure :: tune => Null_tune
        procedure :: try => Null_try
    end type Null_PLMC_Propagator

contains

!implementation Abstract_PLMC_Propagator

    subroutine Abstract_construct(this, generating_algorithms, selector, tuning_period)
        class(Abstract_PLMC_Propagator), intent(out) :: this
        type(Generating_Algorithm_Wrapper), target, intent(in) :: generating_algorithms(:)
        class(Abstract_Tower_Sampler), intent(in) :: selector
        integer, intent(in) :: tuning_period

        this%generating_algorithms => generating_algorithms
        allocate(this%selector, source=selector)
        call check_positive("Abstract_PLMC_Propagator: construct", "tuning_period", tuning_period)
        this%tuning_period = tuning_period
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_PLMC_Propagator), intent(inout) :: this

        call tower_sampler_destroy(this%selector)
        this%generating_algorithms => null()
    end subroutine Abstract_destroy

    subroutine Abstract_reset(this)
        class(Abstract_PLMC_Propagator), intent(inout) :: this

        integer :: nums_algorithms(size(this%generating_algorithms)), i_algorithm

        do i_algorithm = 1, size(this%generating_algorithms)
            call this%generating_algorithms(i_algorithm)%algorithm%reset_selectors()
            nums_algorithms(i_algorithm) = this%generating_algorithms(i_algorithm)%algorithm%&
                get_num_choices()
        end do
        call this%selector%reset(nums_algorithms)
    end subroutine Abstract_reset

    subroutine Abstract_tune(this, i_step)
        class(Abstract_PLMC_Propagator), intent(inout) :: this
        integer, intent(in) :: i_step

        if (i_step /= 0 .and. mod(i_step, this%tuning_period) == 0) then
            call this%reset()
        end if
    end subroutine Abstract_tune

    !> @bug Direct feed in array doesn't work: gfortran bug?
    subroutine Abstract_try(this, observables)
        class(Abstract_PLMC_Propagator), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        integer :: i_choice, i_random

        do i_choice = 1, this%selector%get_num_choices()
            i_random = this%selector%get()
            call this%generating_algorithms(i_random)%algorithm%try(observables)
        end do
    end subroutine Abstract_try

!end implementation Abstract_PLMC_Propagator

!implementation Null_PLMC_Propagator

    subroutine Null_construct(this, generating_algorithms, selector, tuning_period)
        class(Null_PLMC_Propagator), intent(out) :: this
        type(Generating_Algorithm_Wrapper), target, intent(in) :: generating_algorithms(:)
        class(Abstract_Tower_Sampler), intent(in) :: selector
        integer, intent(in) :: tuning_period
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_PLMC_Propagator), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_reset(this)
        class(Null_PLMC_Propagator), intent(inout) :: this
    end subroutine Null_reset

    subroutine Null_tune(this, i_step)
        class(Null_PLMC_Propagator), intent(inout) :: this
        integer, intent(in) :: i_step
    end subroutine Null_tune

    subroutine Null_try(this, observables)
        class(Null_PLMC_Propagator), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables
    end subroutine Null_try

!end implementation Null_PLMC_Propagator

end module classes_plmc_propagator
