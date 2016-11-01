module procedures_selectors_resetters

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_hetero_couples, only: Abstract_Hetero_Couples
use classes_tower_sampler, only: Abstract_Tower_Sampler
use types_component_wrapper, only: Component_Wrapper
use classes_changed_box_size, only: Abstract_Changed_Box_Size
use classes_exchanged_boxes_size, only: Abstract_Exchanged_Boxes_Size

implicit none

private
public :: reset

interface reset
    module procedure :: reset_one_particle, reset_two_particles
    module procedure :: reset_volumes_homo, reset_volumes_hetero
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

    !> @todo What is the best compromise between minval(), maxval(), average() or sum()?
    subroutine reset_two_particles(selectors, couples, components, can_swap)
        class(Abstract_Tower_Sampler), intent(inout) :: selectors(:)
        class(Abstract_Hetero_Couples), intent(in) :: couples(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        logical, intent(in) :: can_swap(:, :)

        integer :: i_box
        integer :: nums_candidates(couples(1)%get_num()), i_couple, ij_couple(2)

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

    subroutine reset_volumes_homo(selectors, changed_boxes_size, components, have_positions)
        class(Abstract_Tower_Sampler), intent(inout) :: selectors(:)
        class(Abstract_Changed_Box_Size), intent(in) :: changed_boxes_size(:)
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
                nums_candidates = ceiling(changed_boxes_size(i_box)%get_frequency_ratio() * &
                    real(nums_candidates, DP))
                nums_candidates = merge(nums_candidates, 1, nums_candidates > 0)
            else
                nums_candidates = 0
            end if

            call selectors(i_box)%reset(nums_candidates)
        end do
    end subroutine reset_volumes_homo

    !> @note The weight of a couple of boxes nums_candidates(i_couple) is the number of particles
    !> inside the 2 boxes.
    !> @todo Test any(have_positions(:, ij_couple)) condition with gfortran
    subroutine reset_volumes_hetero(selector, couples, exchanged_boxes_size, components, &
        have_positions)
        class(Abstract_Tower_Sampler), intent(inout) :: selector
        class(Abstract_Hetero_Couples), intent(in) :: couples
        class(Abstract_Exchanged_Boxes_Size), intent(in) :: exchanged_boxes_size
        type(Component_Wrapper), intent(in) :: components(:, :)
        logical, intent(in) :: have_positions(:, :)

        integer :: nums_candidates(couples%get_num()), i_couple, ij_couple(2)
        integer :: i_component

        do i_couple = 1, size(nums_candidates)
            ij_couple = couples%get(i_couple)
            if (any(have_positions(:, ij_couple(1))) .and. any(have_positions(:, ij_couple(2)))) &
            then
                nums_candidates(i_couple) = 0
                do i_component = 1, size(have_positions, 1)
                    nums_candidates(i_couple) = nums_candidates(i_couple) + &
                        components(i_component, ij_couple(1))%average_num_particles%get() + &
                        components(i_component, ij_couple(2))%average_num_particles%get()
                end do
                nums_candidates(i_couple) = ceiling(exchanged_boxes_size%get_frequency_ratio() * &
                    real(nums_candidates(i_couple), DP))
                nums_candidates(i_couple) = merge(nums_candidates(i_couple), 1, &
                    nums_candidates(i_couple) > 0)
            else
                nums_candidates(i_couple) = 0
            end if
        end do

        call selector%reset(nums_candidates)
    end subroutine reset_volumes_hetero

end module procedures_selectors_resetters
