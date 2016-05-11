module procedures_one_particle_move_factory

use classes_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use classes_one_particle_move, only: Abstract_One_Particle_Move, &
    Concrete_One_Particle_Translation, Concrete_One_Particle_Rotation, Null_One_Particle_Move
use procedures_property_inquirers, only: component_can_translate, component_can_rotate

implicit none

private
public :: create_translation, create_rotation, destroy

contains

    subroutine create_translation(one_particle_translation, physical_model, change_components)
        class(Abstract_One_Particle_Move), allocatable, intent(out) :: one_particle_translation
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Component_Wrapper), intent(in) :: change_components(:)

        class(Abstract_Tower_Sampler), allocatable :: selector_mold
        logical :: can_move(size(change_components))
        integer :: i_component

        do i_component = 1, size(can_move)
            can_move(i_component) = component_can_translate(change_components(i_component)%&
                translated_positions)
        end do

        if (any(can_move)) then
            allocate(Concrete_Tower_Sampler :: selector_mold)
            allocate(Concrete_One_Particle_Translation :: one_particle_translation)
        else
            allocate(Null_Tower_Sampler :: selector_mold)
            allocate(Null_One_Particle_Move :: one_particle_translation)
        end if

        call one_particle_translation%construct(physical_model%environment, physical_model%&
            mixture, physical_model%short_interactions, physical_model%dipolar_interactions, &
            change_components, can_move, selector_mold)
    end subroutine create_translation

    subroutine create_rotation(one_particle_rotation, physical_model, change_components)
        class(Abstract_One_Particle_Move), allocatable, intent(out) :: one_particle_rotation
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
            allocate(Null_One_Particle_Move :: one_particle_rotation)
        end if

        call one_particle_rotation%construct(physical_model%environment, physical_model%mixture, &
            physical_model%short_interactions, physical_model%dipolar_interactions, &
            change_components, can_rotate, selector_mold)
    end subroutine create_rotation

    subroutine destroy(one_particle_move)
        class(Abstract_One_Particle_Move), allocatable, intent(inout) :: one_particle_move

        if (allocated(one_particle_move)) then
            call one_particle_move%destroy()
            deallocate(one_particle_move)
        end if
    end subroutine destroy

end module procedures_one_particle_move_factory
