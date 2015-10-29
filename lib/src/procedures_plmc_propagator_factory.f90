module procedures_plmc_propagator_factory

use class_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithm_Pointer, &
    Metropolis_Algorithms_Wrapper
use class_plmc_propagator, only: Abstract_PLMC_Propagator, Concrete_PLMC_Propagator, &
    Null_PLMC_Propagator

implicit none

private
public :: plmc_propagator_create, plmc_propagator_destroy

interface plmc_propagator_create
    module procedure :: create_propagator
end interface plmc_propagator_create

interface plmc_propagator_destroy
    module procedure :: destroy_propagator
end interface plmc_propagator_destroy

contains

    subroutine create_propagator(propagator, metropolis)
        class(Abstract_PLMC_Propagator), allocatable, intent(out) :: propagator
        type(Metropolis_Algorithms_Wrapper), target, intent(in) :: metropolis

        class(Abstract_Tower_Sampler), allocatable :: selector
        type(Metropolis_Algorithm_Pointer) :: metropolis_algorithms(2)
        integer :: nums_choices(size(metropolis_algorithms))
        integer :: i_choice

        metropolis_algorithms(1)%algorithm => metropolis%one_particle_move
        metropolis_algorithms(2)%algorithm => metropolis%one_particle_rotation
        do i_choice = 1, size(nums_choices)
            nums_choices(i_choice) = metropolis_algorithms(i_choice)%algorithm%get_num_choices()
        end do
        if (sum(nums_choices) == 0) then
            allocate(Null_Tower_Sampler :: selector)
            allocate(Null_PLMC_Propagator :: propagator)
        else
            allocate(Concrete_Tower_Sampler :: selector)
            allocate(Concrete_PLMC_Propagator :: propagator)
        end if
        call selector%construct(nums_choices)
        call propagator%construct(metropolis_algorithms, selector)
        call selector%destroy()
        deallocate(selector)
    end subroutine create_propagator

    subroutine destroy_propagator(propagator)
        class(Abstract_PLMC_Propagator), allocatable, intent(inout) :: propagator

        if (allocated(propagator)) then
            call propagator%destroy()
            deallocate(propagator)
        end if
    end subroutine destroy_propagator

end module procedures_plmc_propagator_factory
