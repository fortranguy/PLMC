module procedures_box_volume_change_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only:tower_sampler_create => create, tower_sampler_destroy =>&
    destroy
use procedures_mixture_properties, only: set_have_positions
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use classes_changed_box_size, only: Abstract_Changed_Box_Size
use classes_generating_algorithm, only: Abstract_Generating_Algorithm, Null_Generating_Algorithm
use classes_box_volume_change, only: Box_Volume_Change
use procedures_environment_inquirers, only: box_size_can_change

implicit none

private
public :: create

contains

    !> @note All particles of a box are considered together. Therefore num_candidates = 1.
    subroutine create(volume_change, physical_model, changed_boxes_size)
        class(Abstract_Generating_Algorithm), allocatable, intent(out) :: volume_change
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        class(Abstract_Changed_Box_Size), intent(in) :: changed_boxes_size(:)

        class(Abstract_Tower_Sampler), allocatable :: selectors(:)
        integer :: num_candidates
        logical :: boxes_size_can_change
        logical :: have_positions(size(physical_model%mixture%components, 1), &
            size(physical_model%mixture%components, 2))

        boxes_size_can_change = all(box_size_can_change(changed_boxes_size))
        call set_have_positions(have_positions, physical_model%mixture%components)
        if (boxes_size_can_change .and. any(have_positions)) then
            allocate(Box_Volume_Change :: volume_change)
        else
            allocate(Null_Generating_Algorithm :: volume_change)
        end if

        num_candidates = 1
        call tower_sampler_create(selectors, size(changed_boxes_size), num_candidates, &
            boxes_size_can_change)
        select type (volume_change)
            type is (Box_Volume_Change)
                call volume_change%construct(physical_model%environment, physical_model%mixture, &
                    physical_model%short_interactions, physical_model%dipolar_interactions_facades,&
                    changed_boxes_size, have_positions, selectors)
            type is (Null_Generating_Algorithm)
            class default
                call error_exit("procedures_box_volume_change_factory: create: volume_change: "//&
                    "unknown type.")
        end select
        call tower_sampler_destroy(selectors)
    end subroutine create

end module procedures_box_volume_change_factory
