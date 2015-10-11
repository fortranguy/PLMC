module class_monte_carlo_propagator

use data_constants, only: num_algorithms
use class_tower_sampler, only: Abstract_Tower_Sampler
use types_metropolis_wrapper, only: Metropolis_Wrapper
implicit none

private

    type, abstract, public :: Abstract_Monte_Carlo_Propagator
    private
        type(Metropolis_Wrapper), pointer :: metropolis
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

    subroutine Abstract_Monte_Carlo_Propagator_construct(this, metropolis, selector)
        class(Abstract_Monte_Carlo_Propagator), intent(out) :: this
        type(Metropolis_Wrapper), target, intent(in) :: metropolis
        class(Abstract_Tower_Sampler), intent(in) :: selector

        integer :: nums_choices(num_algorithms)

        this%metropolis => metropolis
        nums_choices = [metropolis%one_particle_move%get_num_choices(), &
            metropolis%one_particle_rotation%get_num_choices()]
        allocate(this%selector, mold=selector)
        call this%selector%construct(nums_choices)
        this%num_choices = sum(nums_choices)
    end subroutine Abstract_Monte_Carlo_Propagator_construct

    subroutine Abstract_Monte_Carlo_Propagator_destroy(this)
        class(Abstract_Monte_Carlo_Propagator), intent(inout) :: this

        if (allocated(this%selector)) deallocate(this%selector)
        this%metropolis => null()
    end subroutine Abstract_Monte_Carlo_Propagator_destroy

    subroutine Abstract_Monte_Carlo_Propagator_try(this)
        class(Abstract_Monte_Carlo_Propagator), intent(inout) :: this

        integer :: i_choice, i_algorithm

        do i_choice = 1, this%num_choices
            select case (this%selector%get())
                case (1)
                    call this%metropolis%one_particle_move%try()
                case (2)
                    call this%metropolis%one_particle_rotation%try()
            end select
        end do
    end subroutine Abstract_Monte_Carlo_Propagator_try

!implementation Abstract_Monte_Carlo_Propagator

!implementation Null_Monte_Carlo_Propagator

    subroutine Null_Monte_Carlo_Propagator_construct(this, metropolis, selector)
        class(Null_Monte_Carlo_Propagator), intent(out) :: this
        type(Metropolis_Wrapper), target, intent(in) :: metropolis
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
