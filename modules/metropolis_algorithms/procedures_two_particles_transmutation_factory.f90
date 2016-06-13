module procedures_two_particles_transmutation_factory

use classes_hetero_couples, only: Abstract_Hetero_Couples, Null_Hetero_Couples, Full_Hetero_Couples
use classes_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use classes_two_particles_transmutation, only: Abstract_Two_Particles_Transmutation, &
    Concrete_Two_Particles_Transmutation, Null_Two_Particles_Transmutation
use procedures_changes_factory, only: set_can_exchange

implicit none

private
public :: create, destroy

contains

    subroutine create(two_particles_transmutation, physical_model, changes)
        class(Abstract_Two_Particles_Transmutation), allocatable, intent(out) :: &
            two_particles_transmutation
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        class(Abstract_Hetero_Couples), allocatable :: couples
        class(Abstract_Tower_Sampler), allocatable :: selector_mold
        logical :: can_exchange(size(physical_model%mixture%components))

        call set_can_exchange(can_exchange, physical_model%mixture%components)
        if (count(can_exchange) > 1) then
            allocate(Full_Hetero_Couples :: couples)
            allocate(Concrete_Tower_Sampler :: selector_mold)
            allocate(Concrete_Two_Particles_Transmutation :: two_particles_transmutation)
        else
            allocate(Null_Hetero_Couples :: couples)
            allocate(Null_Tower_Sampler :: selector_mold)
            allocate(Null_Two_Particles_Transmutation :: two_particles_transmutation)
        end if

        call couples%construct(size(can_exchange))
        call two_particles_transmutation%construct(physical_model%environment, physical_model%&
            mixture, physical_model%short_interactions, physical_model%dipolar_interactions, &
            changes, can_exchange, couples, selector_mold)
        call couples%destroy()
    end subroutine create

    subroutine destroy(two_particles_transmutation)
        class(Abstract_Two_Particles_Transmutation), allocatable, intent(inout) :: &
            two_particles_transmutation

        if (allocated(two_particles_transmutation)) then
            call two_particles_transmutation%destroy()
            deallocate(two_particles_transmutation)
        end if
    end subroutine destroy

end module procedures_two_particles_transmutation_factory
