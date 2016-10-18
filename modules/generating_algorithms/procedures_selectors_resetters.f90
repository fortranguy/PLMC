module procedures_selectors_resetters

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_hetero_couples, only: Abstract_Hetero_Couples
use classes_tower_sampler, only: Abstract_Tower_Sampler
use types_component_wrapper, only: Component_Wrapper
use types_changed_box_size_wrapper, only: Changed_Box_Size_Line

implicit none

private
public :: reset

interface reset
    module procedure :: reset_one_particle
    module procedure :: reset_two_particles
    module procedure :: reset_volume
end interface reset

contains

    subroutine reset_one_particle(selectors, components, properties)
        class(Abstract_Tower_Sampler), intent(inout) :: selectors(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        logical, intent(in) :: properties(:, :) !! filters

        integer :: i_box
        integer :: nums_candidates(size(properties, 1)), i_component

        do i_box = 1, size(selectors)
            do i_component = 1, size(nums_candidates)
                nums_candidates(i_component) = merge(components(i_component, i_box)%&
                    average_num_particles%get(), 0, properties(i_component, i_box))
            end do
            call selectors(i_box)%reset(nums_candidates)
        end do
    end subroutine reset_one_particle

    !> @todo Multi boxes generalisation
    !> @todo What is the best compromise between minval(), maxval() and average()?
    subroutine reset_two_particles(selectors, couples, components, can_swap)
        class(Abstract_Tower_Sampler), intent(inout) :: selectors(:)
        class(Abstract_Hetero_Couples), intent(in) :: couples(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        logical, intent(in) :: can_swap(:, :)

        integer :: i_box
        integer :: nums_candidates(size(can_swap, 1)), i_couple, ij_couple(2)

        do i_box = 1, size(selectors)
            do i_couple = 1, size(nums_candidates)
                ij_couple = couples(i_box)%get(i_couple)
                nums_candidates(i_couple) = merge(minval([components(ij_couple(1), i_box)%&
                    average_num_particles%get(), &
                    components(ij_couple(2), i_box)%average_num_particles%get()]), 0, &
                    can_swap(ij_couple(1), i_box) .and. can_swap(ij_couple(2), i_box))
            end do
            call selectors(i_box)%reset(nums_candidates)
        end do
    end subroutine reset_two_particles

    !> @note what if volume exchange?
    subroutine reset_volume(selectors, changed_boxes_size, components, have_positions)
        class(Abstract_Tower_Sampler), intent(inout) :: selectors(:)
        type(Changed_Box_Size_Line), intent(in) :: changed_boxes_size(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        logical, intent(in) :: have_positions(:, :)

        integer :: i_box
        integer :: nums_candidates(1), i_component

        do i_box = 1, size(selectors)
            if (any(have_positions(:, i_box))) then
                nums_candidates = 0
                do i_component = 1, size(have_positions, 1)
                    nums_candidates = nums_candidates + components(i_component, i_box)%&
                        average_num_particles%get()
                end do
                nums_candidates = ceiling(changed_boxes_size(i_box)%line(i_box)%changed%&
                    get_frequency_ratio() * real(nums_candidates, DP))
                nums_candidates = merge(nums_candidates, 1, nums_candidates > 0)
            else
                nums_candidates = 0
            end if

            call selectors(i_box)%reset(nums_candidates)
        end do
    end subroutine reset_volume

end module procedures_selectors_resetters
