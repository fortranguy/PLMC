module procedures_plmc_propagation

use procedures_errors, only: error_exit
use class_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use types_metropolis_wrapper, only: num_algorithms, Metropolis_Wrapper
use types_observables_wrapper, only: Mixture_Observables_Wrapper
implicit none

private
public :: plmc_propagator_construct, plmc_propagator_destroy, plmc_propagator_try

    class(Abstract_Tower_Sampler), allocatable :: selector
    integer :: num_choices

contains

    subroutine plmc_propagator_construct(metropolis)
        type(Metropolis_Wrapper), intent(in) :: metropolis

        integer :: nums_choices(num_algorithms)

        nums_choices = [metropolis%one_particle_move%get_num_choices(), &
            metropolis%one_particle_rotation%get_num_choices()]
        num_choices = sum(nums_choices)
        if (num_choices == 0) then
            allocate(Null_Tower_Sampler :: selector)
        else
            allocate(Concrete_Tower_Sampler :: selector)
        end if
        call selector%construct(nums_choices)
    end subroutine plmc_propagator_construct

    subroutine plmc_propagator_destroy()

        if (allocated(selector)) deallocate(selector)
    end subroutine plmc_propagator_destroy

    subroutine plmc_propagator_try(metropolis, observables)
        type(Metropolis_Wrapper), intent(in) :: metropolis
        type(Mixture_Observables_Wrapper), intent(inout) :: observables

        integer :: i_choice

        do i_choice = 1, num_choices
            select case (selector%get())
                case (0)
                case (1)
                    call metropolis%one_particle_move%try(observables)
                case (2)
                    call metropolis%one_particle_rotation%try(observables)
                case default
                    call error_exit("plmc_propagator_try: missing algorithm")
            end select
        end do
    end subroutine plmc_propagator_try

end module procedures_plmc_propagation
