module class_monte_carlo_propagator

use class_tower_sampler, only: Abstract_Tower_Sampler
use class_metropolis_algorithm, only: Abstract_Metropolis_Algorithm, Metropolis_Algorithm_Pointer
implicit none

private

    type, abstract, public :: Abstract_Monte_Carlo_Propagator
    private
        integer :: num_algorithms
        type(Metropolis_Algorithm_Pointer), allocatable :: algorithms(:)
        class(Abstract_Tower_Sampler), allocatable :: selector
        integer :: num_choices
    contains
        procedure :: construct => Abstract_Monte_Carlo_Propagator_construct
        procedure :: destroy => Abstract_Monte_Carlo_Propagator_destroy
        procedure :: add => Abstract_Monte_Carlo_Propagator_add
        procedure :: set => Abstract_Monte_Carlo_Propagator_set
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

    subroutine Abstract_Monte_Carlo_Propagator_construct(this, selector)
        class(Abstract_Monte_Carlo_Propagator), intent(out) :: this
        class(Abstract_Tower_Sampler), intent(in) :: selector

        this%num_algorithms = 0
        allocate(this%selector, mold=selector)
    end subroutine Abstract_Monte_Carlo_Propagator_construct

    subroutine Abstract_Monte_Carlo_Propagator_add(this, algorithm)
        class(Abstract_Monte_Carlo_Propagator), intent(inout) :: this
        class(Abstract_Metropolis_Algorithm), target, intent(in) :: algorithm

        type(Metropolis_Algorithm_Pointer) :: new

        new%ptr => algorithm
        select case (this%num_algorithms)
            case (0)
                this%algorithms = [new]
            case (1:)
                this%algorithms = [this%algorithms, new]
        end select
        this%num_algorithms = this%num_algorithms + 1
    end subroutine Abstract_Monte_Carlo_Propagator_add

    subroutine Abstract_Monte_Carlo_Propagator_set(this)
        class(Abstract_Monte_Carlo_Propagator), intent(inout) :: this

        integer, allocatable :: nums_choices(:)
        integer :: i_algorithm

        nums_choices = [this%algorithms(1)%ptr%get_num_choices()]
        do i_algorithm = 2, this%num_algorithms
            nums_choices = [nums_choices, this%algorithms(i_algorithm)%ptr%get_num_choices()]
        end do
        call this%selector%construct(nums_choices)
        this%num_choices = sum(nums_choices)
        deallocate(nums_choices)
    end subroutine Abstract_Monte_Carlo_Propagator_set

    subroutine Abstract_Monte_Carlo_Propagator_destroy(this)
        class(Abstract_Monte_Carlo_Propagator), intent(inout) :: this

        if (allocated(this%selector)) deallocate(this%selector)
        if (allocated(this%algorithms)) deallocate(this%algorithms)
    end subroutine Abstract_Monte_Carlo_Propagator_destroy

    subroutine Abstract_Monte_Carlo_Propagator_try(this)
        class(Abstract_Monte_Carlo_Propagator), intent(inout) :: this

        integer :: i_choice

        do i_choice = 1, this%num_choices
            call this%algorithms(this%selector%get())%ptr%try()
        end do
    end subroutine Abstract_Monte_Carlo_Propagator_try

!implementation Abstract_Monte_Carlo_Propagator

!implementation Null_Monte_Carlo_Propagator

    subroutine Null_Monte_Carlo_Propagator_construct(this, selector)
        class(Null_Monte_Carlo_Propagator), intent(out) :: this
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
