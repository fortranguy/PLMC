module procedures_two_particles_switch_factory

use classes_hetero_couples, only: Abstract_Hetero_Couples
use procedures_hetero_couples_factory, only: hetero_couples_create_half => create_half, &
    hetero_couples_destroy => destroy
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only:tower_sampler_create => create, tower_sampler_destroy =>&
    destroy
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use classes_two_particles_switch, only: Abstract_Two_Particles_Switch, &
    Concrete_Two_Particles_Switch, Null_Two_Particles_Switch

implicit none

private
public :: create, destroy

contains

    !> @todo update to gemc_components
    subroutine create(two_particles_switch, physical_model)
        class(Abstract_Two_Particles_Switch), allocatable, intent(out) :: two_particles_switch
        type(Physical_Model_Wrapper), intent(in) :: physical_model

        class(Abstract_Hetero_Couples), allocatable :: couples
        class(Abstract_Tower_Sampler), allocatable :: selector
        integer :: num_components

        num_components = size(physical_model%mixture%gemc_components, 1)
        if (num_components > 1) then
            allocate(Concrete_Two_Particles_Switch :: two_particles_switch)
        else
            allocate(Null_Two_Particles_Switch :: two_particles_switch)
        end if

        call hetero_couples_create_half(couples, num_components)
        call tower_sampler_create(selector, couples%get_num(), num_components > 1)
        call two_particles_switch%construct(physical_model%environment, physical_model%mixture%&
            components, physical_model%short_interactions, physical_model%&
            dipolar_interactions_dynamic, physical_model%dipolar_interactions_static, couples, &
            selector)
        call tower_sampler_destroy(selector)
        call hetero_couples_destroy(couples)
    end subroutine create

    subroutine destroy(two_particles_switch)
        class(Abstract_Two_Particles_Switch), allocatable, intent(inout) :: two_particles_switch

        if (allocated(two_particles_switch)) then
            call two_particles_switch%destroy()
            deallocate(two_particles_switch)
        end if
    end subroutine destroy

end module procedures_two_particles_switch_factory
