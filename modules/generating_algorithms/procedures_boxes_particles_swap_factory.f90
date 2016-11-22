module procedures_boxes_particles_swap_factory

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
use classes_boxes_particles_swap, only: Boxes_Particles_Swap

implicit none

private
public :: create

contains

    subroutine create(particles_swap, physical_model, changes)
        class(Abstract_Generating_Algorithm), allocatable, intent(out) :: particles_swap
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes

        class(Abstract_Hetero_Couples), allocatable :: box_couples, component_couples(:)
        class(Abstract_Tower_Sampler), allocatable :: boxes_selector, components_selectors(:)
        integer, allocatable :: nums_box_couples(:)
        logical :: can_translate(size(changes%components, 1), size(changes%components, 2))

        call set_can_translate(can_translate, changes%components)
        if (size(can_translate, 2) > 1 .and. any(can_translate)) then
            allocate(Boxes_Particles_Swap :: particles_swap)
        else
            allocate(Null_Generating_Algorithm :: particles_swap)
        end if

        call hetero_couples_create_full(box_couples, size(can_translate, 2))
        call tower_sampler_create(boxes_selector, box_couples%get_num(), size(can_translate, 2) > 1)
        allocate(nums_box_couples(box_couples%get_num()))
        nums_box_couples = 1
        call boxes_selector%reset(nums_box_couples)
        call hetero_couples_create_full(component_couples, box_couples%get_num(), &
            size(can_translate, 1))
        call tower_sampler_create(components_selectors, box_couples%get_num(), &
            component_couples(1)%get_num(), any(can_translate))

        select type (particles_swap)
            type is (Boxes_Particles_Swap)
                call particles_swap%construct(physical_model%environment, physical_model%mixture, &
                    physical_model%short_interactions, physical_model%dipolar_interactions_dynamic,&
                    physical_model%dipolar_interactions_static, can_translate, box_couples, &
                    component_couples, boxes_selector, components_selectors)
            type is (Null_Generating_Algorithm)
            class default
                call error_exit("procedures_boxes_particles_swap_factory: create: particles_swap: "&
                    //"unknown type.")
        end select

        call tower_sampler_destroy(components_selectors)
        call hetero_couples_destroy(component_couples)
        call tower_sampler_destroy(boxes_selector)
        call hetero_couples_destroy(box_couples)
    end subroutine create

end module procedures_boxes_particles_swap_factory
