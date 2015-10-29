module procedures_plmc_propagator

use class_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithm_Pointer, &
    Metropolis_Algorithms_Wrapper
use types_observables_wrapper, only: Observables_Wrapper

implicit none

private
public :: plmc_propagator_create, plmc_propagator_destroy, plmc_propagator_try

    type(Metropolis_Algorithm_Pointer) :: algorithms(2)
    class(Abstract_Tower_Sampler), allocatable :: selector

contains

    subroutine plmc_propagator_create(metropolis_algorithms)
        type(Metropolis_Algorithms_Wrapper), target, intent(in) :: metropolis_algorithms

        integer :: nums_choices(size(algorithms))
        integer :: i_choice

        algorithms(1)%algorithm => metropolis_algorithms%one_particle_move
        algorithms(2)%algorithm => metropolis_algorithms%one_particle_rotation
        do i_choice = 1, size(nums_choices)
            nums_choices(i_choice) = algorithms(i_choice)%algorithm%get_num_choices()
        end do
        if (sum(nums_choices) == 0) then
            allocate(Null_Tower_Sampler :: selector)
        else
            allocate(Concrete_Tower_Sampler :: selector)
        end if
        call selector%construct(nums_choices)
    end subroutine plmc_propagator_create

    subroutine plmc_propagator_destroy()

        integer :: i_algorithm

        if (allocated(selector)) then
            call selector%destroy()
            deallocate(selector)
        end if
        do i_algorithm = size(algorithms), 1, -1
            algorithms(i_algorithm)%algorithm => null()
        end do
    end subroutine plmc_propagator_destroy

    subroutine plmc_propagator_try(observables)
        type(Observables_Wrapper), intent(inout) :: observables

        integer :: i_choice

        do i_choice = 1, selector%get_num_choices()
            call algorithms(selector%get())%algorithm%try(observables)
        end do
    end subroutine plmc_propagator_try

end module procedures_plmc_propagator
