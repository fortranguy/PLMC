module classes_plmc_propagator

use procedures_checks, only: check_positive
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only: tower_sampler_destroy => destroy
use types_component_wrapper, only: Component_Wrapper
use types_generating_algorithms_wrapper, only: Generating_Algorithm_Pointer, &
    Generating_Algorithms_Wrapper
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper

implicit none

private

    type, abstract, public :: Abstract_PLMC_Propagator
    private
        type(Component_Wrapper), pointer :: components(:, :) => null()
        type(Generating_Algorithm_Pointer), allocatable :: generating_algorithms(:)
        integer :: accumulation_period
        class(Abstract_Tower_Sampler), allocatable :: selector
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset_selector => Abstract_reset_selector
        procedure :: tune => Abstract_tune
        procedure :: try => Abstract_try
        procedure, private :: set_accumulation_period => Abstract_set_accumulation_period
    end type Abstract_PLMC_Propagator

    type, extends(Abstract_PLMC_Propagator), public :: Concrete_PLMC_Propagator

    end type Concrete_PLMC_Propagator

    type, extends(Abstract_PLMC_Propagator), public :: Null_PLMC_Propagator
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: reset_selector => Null_reset_selector
        procedure :: tune => Null_tune
        procedure :: try => Null_try
    end type Null_PLMC_Propagator

contains

!implementation Abstract_PLMC_Propagator

    subroutine Abstract_construct(this, components, generating_algorithms, selector)
        class(Abstract_PLMC_Propagator), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: components(:, :)
        type(Generating_Algorithm_Pointer), target, intent(in) :: generating_algorithms(:)
        class(Abstract_Tower_Sampler), intent(in) :: selector

        this%components => components
        call this%set_accumulation_period()
        allocate(this%generating_algorithms, source=generating_algorithms)
        allocate(this%selector, source=selector)
    end subroutine Abstract_construct

    !> @note The choice of maxval is arbitrary. It looks like a dirty hack.
    subroutine Abstract_set_accumulation_period(this)
        class(Abstract_PLMC_Propagator), intent(inout) :: this

        integer :: accumulation_periods(size(this%components, 1), size(this%components, 2))
        integer :: i_box, i_component

        do i_box = 1, size(this%components, 2)
            do i_component = 1, size(this%components, 1)
                accumulation_periods(i_component, i_box) = this%components(i_component, i_box)%&
                    average_num_particles%get_accumulation_period()
            end do
        end do
        call check_positive("Abstract_PLMC_Propagator: set_accumulation_period", &
            "maxval(accumulation_periods)", maxval(accumulation_periods))
        this%accumulation_period = maxval(accumulation_periods)
    end subroutine Abstract_set_accumulation_period

    subroutine Abstract_destroy(this)
        class(Abstract_PLMC_Propagator), intent(inout) :: this

        call tower_sampler_destroy(this%selector)
        if (allocated(this%generating_algorithms)) then
            deallocate(this%generating_algorithms)
        end if
        this%components => null()
    end subroutine Abstract_destroy

    subroutine Abstract_reset_selector(this)
        class(Abstract_PLMC_Propagator), intent(inout) :: this

        integer :: nums_algorithms(size(this%generating_algorithms)), i_algorithm

        do i_algorithm = 1, size(this%generating_algorithms)
            call this%generating_algorithms(i_algorithm)%algorithm%reset_selector()
            nums_algorithms(i_algorithm) = this%generating_algorithms(i_algorithm)%algorithm%&
                get_num_choices()
        end do
        call this%selector%reset(nums_algorithms)
    end subroutine Abstract_reset_selector

    subroutine Abstract_tune(this, i_step)
        class(Abstract_PLMC_Propagator), intent(inout) :: this
        integer, intent(in) :: i_step

        integer :: i_box, i_component

        do i_box = 1, size(this%components, 2)
            do i_component = 1, size(this%components, 1)
                call this%components(i_component, i_box)%average_num_particles%accumulate(i_step)
            end do
        end do

        if (mod(i_step, this%accumulation_period) == 0) then
            call this%reset_selector()
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

    subroutine Null_construct(this, components, generating_algorithms, selector)
        class(Null_PLMC_Propagator), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: components(:, :)
        type(Generating_Algorithm_Pointer), target, intent(in) :: generating_algorithms(:)
        class(Abstract_Tower_Sampler), intent(in) :: selector
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_PLMC_Propagator), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_reset_selector(this)
        class(Null_PLMC_Propagator), intent(inout) :: this
    end subroutine Null_reset_selector

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
