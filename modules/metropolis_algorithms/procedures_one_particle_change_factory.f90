module procedures_one_particle_change_factory

use classes_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use classes_one_particle_change, only: Abstract_One_Particle_Change, &
    Concrete_One_Particle_Move, Concrete_One_Particle_Rotation, Null_One_Particle_Change
use procedures_property_inquirers, only: component_can_move, component_can_rotate

implicit none

private
public :: create_move, create_rotation, destroy

contains

    subroutine create_move(one_particle_move, physical_model, change_components)
        class(Abstract_One_Particle_Change), allocatable, intent(out) :: one_particle_move
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Component_Wrapper), intent(in) :: change_components(:)

        class(Abstract_Tower_Sampler), allocatable :: selector_mold
        logical :: can_move(size(change_components))
        integer :: i_component

        do i_component = 1, size(can_move)
            can_move(i_component) = component_can_move(change_components(i_component)%&
                moved_positions)
        end do

        if (any(can_move)) then
            allocate(Concrete_Tower_Sampler :: selector_mold)
            allocate(Concrete_One_Particle_Move :: one_particle_move)
        else
            allocate(Null_Tower_Sampler :: selector_mold)
            allocate(Null_One_Particle_Change :: one_particle_move)
        end if

        call one_particle_move%construct(physical_model%environment, physical_model%mixture, &
            physical_model%short_interactions, physical_model%dipolar_interactions, &
            change_components, can_move, selector_mold)
    end subroutine create_move

    subroutine create_rotation(one_particle_rotation, physical_model, change_components)
        class(Abstract_One_Particle_Change), allocatable, intent(out) :: one_particle_rotation
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Component_Wrapper), intent(in) :: change_components(:)

        class(Abstract_Tower_Sampler), allocatable :: selector_mold
        logical :: can_rotate(size(change_components))
        integer :: i_component

        do i_component = 1, size(can_rotate)
            can_rotate(i_component) = component_can_rotate(change_components(i_component)%&
                rotated_orientations)
        end do

        if (any(can_rotate)) then
            allocate(Concrete_Tower_Sampler :: selector_mold)
            allocate(Concrete_One_Particle_Rotation :: one_particle_rotation)
        else
            allocate(Null_Tower_Sampler :: selector_mold)
            allocate(Null_One_Particle_Change :: one_particle_rotation)
        end if

        call one_particle_rotation%construct(physical_model%environment, physical_model%mixture, &
            physical_model%short_interactions, physical_model%dipolar_interactions, &
            change_components, can_rotate, selector_mold)
    end subroutine create_rotation

    subroutine destroy(one_particle_change)
        class(Abstract_One_Particle_Change), allocatable, intent(inout) :: one_particle_change

        if (allocated(one_particle_change)) then
            call one_particle_change%destroy()
            deallocate(one_particle_change)
        end if
    end subroutine destroy

end module procedures_one_particle_change_factory
