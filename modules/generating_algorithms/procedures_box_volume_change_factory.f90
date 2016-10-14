module procedures_box_volume_change_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only:tower_sampler_create => create, tower_sampler_destroy =>&
    destroy
use procedures_mixture_factory, only: set_have_positions
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use classes_changed_box_size, only: Abstract_Changed_Box_Size
use classes_box_volume_change, only: Abstract_Box_Volume_Change, Concrete_Box_Volume_Change, &
    Null_Box_Volume_Change
use procedures_environment_inquirers, only: property_box_size_can_change => box_size_can_change

implicit none

private
public :: create, destroy

contains

    !> @note All particles of a box are considered together. Therefore num_candidates = 1.
    subroutine create(box_volume_change, physical_model, changed_box_size)
        class(Abstract_Box_Volume_Change), allocatable, intent(out) :: box_volume_change
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        class(Abstract_Changed_Box_Size), intent(in) :: changed_box_size

        class(Abstract_Tower_Sampler), allocatable :: selector
        integer :: num_candidates
        logical :: box_size_can_change
        logical :: have_positions(size(physical_model%mixture%gemc_components, 1))

        box_size_can_change = property_box_size_can_change(changed_box_size)
        call set_have_positions(have_positions, physical_model%mixture%components)
        if (box_size_can_change) then
            allocate(Concrete_Box_Volume_Change :: box_volume_change)
        else
            allocate(Null_Box_Volume_Change :: box_volume_change)
        end if

        num_candidates = 1
        call tower_sampler_create(selector, num_candidates, box_size_can_change)
        call box_volume_change%construct(physical_model%environment, physical_model%mixture, &
            physical_model%short_interactions, physical_model%dipolar_interactions_facade, &
            changed_box_size, have_positions, selector)
        call tower_sampler_destroy(selector)
    end subroutine create

    subroutine destroy(box_volume_change)
        class(Abstract_Box_Volume_Change), allocatable, intent(inout) :: box_volume_change

        if (allocated(box_volume_change)) then
            call box_volume_change%destroy()
            deallocate(box_volume_change)
        end if
    end subroutine destroy

end module procedures_box_volume_change_factory
