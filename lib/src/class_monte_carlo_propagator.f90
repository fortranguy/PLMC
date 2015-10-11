module class_monte_carlo_propagator

use class_tower_sampler, only: Abstract_Tower_Sampler
use class_metropolis_algorithm, only: Abstract_Metropolis_Algorithm
implicit none

private

    type, abstract, public :: Abstract_Monte_Carlo_Propagator
    private
        class(Abstract_Metropolis_Algorithm), pointer :: algorithms(:)
        class(Abstract_Tower_Sampler), allocatable :: selector
        integer :: num_choices
    contains
        procedure :: construct => Abstract_Monte_Carlo_Propagator_construct
        procedure :: destroy => Abstract_Monte_Carlo_Propagator_destroy
        procedure :: try => Abstract_Monte_Carlo_Propagator_try
    end type Abstract_Monte_Carlo_Propagator

    type, extends(Abstract_Monte_Carlo_Propagator), public :: Concrete_Monte_Carlo_Propagator

    end type Concrete_Monte_Carlo_Propagator

    type, extends(Abstract_Monte_Carlo_Propagator), public :: Null_Monte_Carlo_Propagator
    contains
        procedure :: construct => Null_Monte_Carlo_Propagator_construct
        procedure :: destroy => Null_Monte_Carlo_Propagator_destroy
        procedure :: try => Null_Monte_Carlo_Propagator_try
    end type Null_Monte_Carlo_Propagator

contains

!implementation Abstract_Monte_Carlo_Propagator

    subroutine Abstract_Monte_Carlo_Propagator_construct(this, algorithms, selector)
        class(Abstract_Monte_Carlo_Propagator), intent(out) :: this
        class(Abstract_Metropolis_Algorithm), target, intent(in) :: algorithms(:)
        class(Abstract_Tower_Sampler), intent(in) :: selector

        integer, allocatable :: nums_choices(:)
        integer :: i_algorithm

        this%algorithms => algorithms
        allocate(this%selector, mold=selector)
        write(*, *) "algorithms(1)%get_num_choices()", this%algorithms(1)%get_num_choices()
        nums_choices = [this%algorithms(1)%get_num_choices()]
        do i_algorithm = 2, size(this%algorithms)
            nums_choices = [nums_choices, this%algorithms(i_algorithm)%get_num_choices()]
        end do
        call this%selector%construct(nums_choices)
        this%num_choices = sum(nums_choices)
        deallocate(nums_choices)
    end subroutine Abstract_Monte_Carlo_Propagator_construct

    subroutine Abstract_Monte_Carlo_Propagator_destroy(this)
        class(Abstract_Monte_Carlo_Propagator), intent(inout) :: this

        if (allocated(this%selector)) deallocate(this%selector)
        this%algorithms => null()
    end subroutine Abstract_Monte_Carlo_Propagator_destroy

    subroutine Abstract_Monte_Carlo_Propagator_try(this)
        class(Abstract_Monte_Carlo_Propagator), intent(inout) :: this

        integer :: i_choice

        do i_choice = 1, this%num_choices
            call this%algorithms(this%selector%get())%try()
        end do
    end subroutine Abstract_Monte_Carlo_Propagator_try

!implementation Abstract_Monte_Carlo_Propagator

!implementation Null_Monte_Carlo_Propagator

    subroutine Null_Monte_Carlo_Propagator_construct(this, algorithms, selector)
        class(Null_Monte_Carlo_Propagator), intent(out) :: this
        class(Abstract_Metropolis_Algorithm), target, intent(in) :: algorithms(:)
        class(Abstract_Tower_Sampler), intent(in) :: selector
    end subroutine Null_Monte_Carlo_Propagator_construct

    subroutine Null_Monte_Carlo_Propagator_destroy(this)
        class(Null_Monte_Carlo_Propagator), intent(inout) :: this
    end subroutine Null_Monte_Carlo_Propagator_destroy

    subroutine Null_Monte_Carlo_Propagator_try(this)
        class(Null_Monte_Carlo_Propagator), intent(inout) :: this
    end subroutine Null_Monte_Carlo_Propagator_try

!implementation Null_Monte_Carlo_Propagator

end module class_monte_carlo_propagator
