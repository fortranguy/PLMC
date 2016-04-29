module procedures_two_particles_switch_factory

use classes_hetero_couples, only: Abstract_Hetero_Couples, Null_Hetero_Couples, &
    Concrete_Hetero_Couples
use classes_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use classes_two_particles_switch, only: Abstract_Two_Particles_Switch, &
    Concrete_Two_Particles_Switch, Null_Two_Particles_Switch

implicit none

private
public :: create, destroy, set

contains

    subroutine create(two_particles_switch, environment, components, short_interactions, &
        dipolar_interactions)
        class(Abstract_Two_Particles_Switch), allocatable, intent(out) :: two_particles_switch
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        if (size(components) > 1) then
            allocate(Concrete_Two_Particles_Switch :: two_particles_switch)
        else
            allocate(Null_Two_Particles_Switch :: two_particles_switch)
        end if

        call two_particles_switch%construct(environment, components, short_interactions, &
            dipolar_interactions)
    end subroutine create

    subroutine destroy(two_particles_switch)
        class(Abstract_Two_Particles_Switch), allocatable, intent(inout) :: two_particles_switch

        if (allocated(two_particles_switch)) then
            call two_particles_switch%destroy()
            deallocate(two_particles_switch)
        end if
    end subroutine destroy

    subroutine set(two_particles_switch, components)
        class(Abstract_Two_Particles_Switch), intent(inout) :: two_particles_switch
        type(Component_Wrapper), intent(in) :: components(:)

        class(Abstract_Hetero_Couples), allocatable :: couples
        class(Abstract_Tower_Sampler), allocatable :: selector

        integer, allocatable :: nums_candidates(:)
        integer :: i_candidate, ij_couple(2)

        if (size(components) > 1) then
            allocate(Concrete_Hetero_Couples :: couples)
            allocate(Concrete_Tower_Sampler :: selector)
        else
            allocate(Null_Hetero_Couples :: couples)
            allocate(Null_Tower_Sampler :: selector)
        end if
        call couples%construct(size(components))

        allocate(nums_candidates(couples%get_num_indices()))
        do i_candidate = 1, size(nums_candidates)
            ij_couple = couples%get(i_candidate)
            nums_candidates(i_candidate) = minval([components(ij_couple(1))%average_number%get(), &
                components(ij_couple(2))%average_number%get()])
                !What is the best compromise: minval, maxval, average?
        end do
        call selector%construct(nums_candidates)

        call two_particles_switch%set_couples_and_selector(couples, selector)
        call selector%destroy()
        call couples%destroy()
    end subroutine set

end module procedures_two_particles_switch_factory
