module procedures_one_particle_move_factory

use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only:tower_sampler_create => create, tower_sampler_destroy =>&
    destroy
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use procedures_changes_factory, only: set_can_translate, set_can_rotate
use classes_one_particle_move, only: Abstract_One_Particle_Move, &
    Concrete_One_Particle_Translation, Concrete_One_Particle_Rotation, Null_One_Particle_Move

implicit none

private
public :: create_translation, create_rotation, destroy

contains

    subroutine create_translation(one_particle_translation, physical_model, changes_components)
        class(Abstract_One_Particle_Move), allocatable, intent(out) :: one_particle_translation
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Component_Wrapper), intent(in) :: changes_components(:)

        class(Abstract_Tower_Sampler), allocatable :: selector
        logical :: can_translate(size(changes_components))

        call set_can_translate(can_translate, changes_components)
        if (any(can_translate)) then
            allocate(Concrete_One_Particle_Translation :: one_particle_translation)
        else
            allocate(Null_One_Particle_Move :: one_particle_translation)
        end if

        call tower_sampler_create(selector, size(can_translate), any(can_translate))
        call one_particle_translation%construct(physical_model%environment, physical_model%&
            mixture, physical_model%short_interactions, physical_model%&
            dipolar_interactions_dynamic, physical_model%dipolar_interactions_static, &
            changes_components, can_translate, selector)
        call tower_sampler_destroy(selector)
    end subroutine create_translation

    subroutine create_rotation(one_particle_rotation, physical_model, changes_components)
        class(Abstract_One_Particle_Move), allocatable, intent(out) :: one_particle_rotation
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Component_Wrapper), intent(in) :: changes_components(:)

        class(Abstract_Tower_Sampler), allocatable :: selector
        logical :: can_rotate(size(changes_components))

        call set_can_rotate(can_rotate, changes_components)
        if (any(can_rotate)) then
            allocate(Concrete_One_Particle_Rotation :: one_particle_rotation)
        else
            allocate(Null_One_Particle_Move :: one_particle_rotation)
        end if

        call tower_sampler_create(selector, size(can_rotate), any(can_rotate))
        call one_particle_rotation%construct(physical_model%environment, physical_model%mixture, &
            physical_model%short_interactions, physical_model%dipolar_interactions_dynamic, &
            physical_model%dipolar_interactions_static, changes_components, can_rotate, &
            selector)
        call tower_sampler_destroy(selector)
    end subroutine create_rotation

    subroutine destroy(one_particle_move)
        class(Abstract_One_Particle_Move), allocatable, intent(inout) :: one_particle_move

        if (allocated(one_particle_move)) then
            call one_particle_move%destroy()
            deallocate(one_particle_move)
        end if
    end subroutine destroy

end module procedures_one_particle_move_factory
