module procedures_box_particle_move_factory

use procedures_errors, only: error_exit
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only:tower_sampler_create => create, tower_sampler_destroy =>&
    destroy
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use procedures_changes_factory, only: set_can_translate, set_can_rotate
use classes_generating_algorithm, only: Abstract_Generating_Algorithm, Null_Generating_Algorithm
use classes_box_particle_move, only: Box_Particle_Translation, Box_Particle_Rotation

implicit none

private
public :: create_translation, create_rotation

contains

    subroutine create_translation(particle_translation, physical_model, changes_components)
        class(Abstract_Generating_Algorithm), allocatable, intent(out) :: particle_translation
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Component_Wrapper), intent(in) :: changes_components(:, :)

        class(Abstract_Tower_Sampler), allocatable :: selectors(:)
        logical :: can_translate(size(changes_components, 1), size(changes_components, 2))

        call set_can_translate(can_translate, changes_components)
        if (any(can_translate)) then
            allocate(Box_Particle_Translation :: particle_translation)
        else
            allocate(Null_Generating_Algorithm :: particle_translation)
        end if

        call tower_sampler_create(selectors, size(can_translate, 2), size(can_translate, 1), &
            any(can_translate))
        select type (particle_translation)
            type is (Box_Particle_Translation)
                call particle_translation%construct(physical_model%environment, physical_model%&
                    mixture, physical_model%short_interactions, physical_model%&
                    dipolar_interactions_dynamic, physical_model%dipolar_interactions_static, &
                    changes_components, can_translate, selectors)
            type is (Null_Generating_Algorithm)
            class default
                call error_exit("procedures_box_particle_move_factory: create_translation: "//&
                    "particle_translation: unknown type.")
        end select
        call tower_sampler_destroy(selectors)
    end subroutine create_translation

    subroutine create_rotation(particle_rotation, physical_model, changes_components)
        class(Abstract_Generating_Algorithm), allocatable, intent(out) :: particle_rotation
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Component_Wrapper), intent(in) :: changes_components(:, :)

        class(Abstract_Tower_Sampler), allocatable :: selectors(:)
        logical :: can_rotate(size(changes_components, 1), size(changes_components, 2))

        call set_can_rotate(can_rotate, changes_components)
        if (any(can_rotate)) then
            allocate(Box_Particle_Rotation :: particle_rotation)
        else
            allocate(Null_Generating_Algorithm :: particle_rotation)
        end if

        call tower_sampler_create(selectors, size(can_rotate, 2), size(can_rotate, 1), &
            any(can_rotate))
        select type (particle_rotation)
            type is (Box_Particle_Rotation)
                call particle_rotation%construct(physical_model%environment, &
                    physical_model%mixture, physical_model%short_interactions, &
                    physical_model%dipolar_interactions_dynamic, &
                    physical_model%dipolar_interactions_static, changes_components, &
                    can_rotate, selectors)
            type is (Null_Generating_Algorithm)
            class default
                call error_exit("procedures_box_particle_move_factory: create_rotation: "//&
                    "particle_rotation: unknown type.")
        end select
        call tower_sampler_destroy(selectors)
    end subroutine create_rotation

end module procedures_box_particle_move_factory
