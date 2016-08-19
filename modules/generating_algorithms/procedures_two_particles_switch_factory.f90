module procedures_two_particles_switch_factory

use classes_hetero_couples, only: Abstract_Hetero_Couples, Null_Hetero_Couples, Half_Hetero_Couples
use classes_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use classes_two_particles_switch, only: Abstract_Two_Particles_Switch, &
    Concrete_Two_Particles_Switch, Null_Two_Particles_Switch

implicit none

private
public :: create, destroy

contains

    subroutine create(two_particles_switch, physical_model)
        class(Abstract_Two_Particles_Switch), allocatable, intent(out) :: two_particles_switch
        type(Physical_Model_Wrapper), intent(in) :: physical_model

        class(Abstract_Hetero_Couples), allocatable :: couples
        class(Abstract_Tower_Sampler), allocatable :: selector_mold

        if (size(physical_model%mixture%components) > 1) then
            allocate(Half_Hetero_Couples :: couples)
            allocate(Concrete_Tower_Sampler :: selector_mold)
            allocate(Concrete_Two_Particles_Switch :: two_particles_switch)
        else
            allocate(Null_Hetero_Couples :: couples)
            allocate(Null_Tower_Sampler :: selector_mold)
            allocate(Null_Two_Particles_Switch :: two_particles_switch)
        end if

        call couples%construct(size(physical_model%mixture%components))
        call two_particles_switch%construct(physical_model%environment, physical_model%mixture%&
            components, physical_model%short_interactions, physical_model%dipolar_interactions, &
            couples, selector_mold)
        call couples%destroy()
    end subroutine create

    subroutine destroy(two_particles_switch)
        class(Abstract_Two_Particles_Switch), allocatable, intent(inout) :: two_particles_switch

        if (allocated(two_particles_switch)) then
            call two_particles_switch%destroy()
            deallocate(two_particles_switch)
        end if
    end subroutine destroy

end module procedures_two_particles_switch_factory
