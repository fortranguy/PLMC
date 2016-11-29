module procedures_boxes_particle_teleportation_factory

use procedures_errors, only: error_exit
use classes_hetero_couples, only: Abstract_Hetero_Couples
use procedures_hetero_couples_factory, only: hetero_couples_create_full => create_full, &
    hetero_couples_destroy => destroy
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only:tower_sampler_create => create, tower_sampler_destroy =>&
    destroy
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use procedures_changes_properties, only: set_can_translate
use classes_generating_algorithm, only: Abstract_Generating_Algorithm, Null_Generating_Algorithm
use classes_boxes_particle_teleportation, only: Boxes_Particle_Teleportation

implicit none

private
public :: create

contains

    !> @todo better condition?
    subroutine create(particle_teleportation, physical_model, changes)
        class(Abstract_Generating_Algorithm), allocatable, intent(out) :: particle_teleportation
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        class(Abstract_Hetero_Couples), allocatable :: box_couples
        class(Abstract_Tower_Sampler), allocatable :: boxes_selector, component_selectors(:)
        integer, allocatable :: nums_box_couples(:)
        logical :: can_translate(size(changes%components, 1), size(changes%components, 2))

        call set_can_translate(can_translate, changes%components)
        if (size(can_translate, 2) > 1 .and. any(can_translate)) then
            allocate(Boxes_Particle_Teleportation :: particle_teleportation)
        else
            allocate(Null_Generating_Algorithm :: particle_teleportation)
        end if

        call hetero_couples_create_full(box_couples, size(can_translate, 2))
        call tower_sampler_create(boxes_selector, box_couples%get_num(), size(can_translate, 2) > 1)
        allocate(nums_box_couples(box_couples%get_num()))
        nums_box_couples = 1
        call boxes_selector%reset(nums_box_couples)
        call tower_sampler_create(component_selectors, size(can_translate, 2), &
            size(can_translate, 1), any(can_translate))

        select type (particle_teleportation)
            type is (Boxes_Particle_Teleportation)
                call particle_teleportation%construct(physical_model%environment, physical_model%&
                    mixture, physical_model%short_interactions, physical_model%&
                    dipolar_interactions_dynamic, physical_model%dipolar_interactions_static, &
                    changes%random_positions, can_translate, box_couples, boxes_selector, &
                    component_selectors)
            type is (Null_Generating_Algorithm)
            class default
                call error_exit("procedures_boxes_particle_teleportation_factory: create: "//&
                    "particle_teleportation: unknown type.")
        end select

        call tower_sampler_destroy(component_selectors)
        call tower_sampler_destroy(boxes_selector)
        call hetero_couples_destroy(box_couples)
    end subroutine create

end module procedures_boxes_particle_teleportation_factory
