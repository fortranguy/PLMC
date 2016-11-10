module procedures_selectors_resetters

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_hetero_couples, only: Abstract_Hetero_Couples
use classes_tower_sampler, only: Abstract_Tower_Sampler
use types_component_wrapper, only: Component_Wrapper
use classes_changed_box_size, only: Abstract_Changed_Box_Size
use classes_exchanged_boxes_size, only: Exchanged_Boxes_Size_Line

implicit none

private
public :: reset

interface reset
    module procedure :: reset_box_particle, reset_box_particles_swap, reset_boxes_particles_swap
    module procedure :: reset_box_volume_change, reset_boxes_volume_exchange
end interface reset

contains

    subroutine reset_box_particle(selectors, components, properties)
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
    end subroutine reset_box_particle

    !> @todo What is the best compromise between minval(), maxval(), average() or sum()?
    subroutine reset_box_particles_swap(selectors, couples, components, can_swap)
        class(Abstract_Tower_Sampler), intent(inout) :: selectors(:)
        class(Abstract_Hetero_Couples), intent(in) :: couples(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        logical, intent(in) :: can_swap(:, :)

        integer :: i_box
        integer :: i_couple, ij_components(2)
        integer :: nums_candidates(couples(1)%get_num()), nums_particles(size(ij_components))

        do i_box = 1, size(selectors)
            do i_couple = 1, size(nums_candidates)
                ij_components = couples(i_box)%get(i_couple)
                nums_particles = [components(ij_components(1), i_box)%average_num_particles%get(), &
                    components(ij_components(2), i_box)%average_num_particles%get()]
                nums_candidates(i_couple) = merge(maxval(nums_particles), 0, &
                        can_swap(ij_components(1), i_box) .and. can_swap(ij_components(2), i_box))
            end do
            call selectors(i_box)%reset(nums_candidates)
        end do
    end subroutine reset_box_particles_swap

    subroutine reset_boxes_particles_swap(selectors, box_couples, component_couples, components, &
        can_translate)
        class(Abstract_Tower_Sampler), intent(inout) :: selectors(:)
        class(Abstract_Hetero_Couples), intent(in) :: box_couples, component_couples(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        logical, intent(in) :: can_translate(:, :)

        integer :: box_i_couple, ij_boxes(2)
        integer :: component_i_couple, ij_components(2), i_component
        integer :: boxes_nums_particles(size(components, 1)), num_particles
        integer :: nums_candidates(component_couples(1)%get_num())

        do box_i_couple = 1, box_couples%get_num()
            ij_boxes = box_couples%get(box_i_couple)
            do i_component = 1, size(boxes_nums_particles)
                boxes_nums_particles(i_component) = &
                    components(i_component, ij_boxes(1))%average_num_particles%get() + &
                    components(i_component, ij_boxes(2))%average_num_particles%get()
            end do
            do component_i_couple = 1, component_couples(box_i_couple)%get_num()
                ij_components = component_couples(box_i_couple)%get(component_i_couple)
                num_particles = boxes_nums_particles(ij_components(1)) + &
                    boxes_nums_particles(ij_components(2))
                nums_candidates(component_i_couple) = merge(num_particles/2, 0, &
                    can_translate(ij_components(1), ij_boxes(1)) .and. &
                    can_translate(ij_components(2), ij_boxes(2)))
            end do
            call selectors(box_i_couple)%reset(nums_candidates)
        end do
    end subroutine reset_boxes_particles_swap

    subroutine reset_box_volume_change(selectors, changed_boxes_size, components, have_positions)
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
    end subroutine reset_box_volume_change

    !> @note The weight of a couple of boxes nums_candidates(i_couple) is the number of particles
    !> inside the 2 boxes.
    !> @todo Test any(have_positions(:, ij_boxes)) condition with gfortran
    !> @note direct feed of maxval / minval causes a bug.
    !> @todo Send a bug report to gfortran
    subroutine reset_boxes_volume_exchange(selector, couples, exchanged_boxes_size, components, &
        have_positions)
        class(Abstract_Tower_Sampler), intent(inout) :: selector
        class(Abstract_Hetero_Couples), intent(in) :: couples
        type(Exchanged_Boxes_Size_Line), intent(in) :: exchanged_boxes_size(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        logical, intent(in) :: have_positions(:, :)

        integer :: nums_candidates(couples%get_num()), i_couple, ij_boxes(2), i_box, j_box
        integer :: i_component

        do i_couple = 1, size(nums_candidates)
            ij_boxes = couples%get(i_couple)
            j_box = maxval(ij_boxes)
            i_box = minval(ij_boxes)
            if (any(have_positions(:, i_box)) .and. any(have_positions(:, j_box))) then
                nums_candidates(i_couple) = 0
                do i_component = 1, size(have_positions, 1)
                    nums_candidates(i_couple) = nums_candidates(i_couple) + &
                        components(i_component, ij_boxes(1))%average_num_particles%get() + &
                        components(i_component, ij_boxes(2))%average_num_particles%get()
                end do
                nums_candidates(i_couple) = ceiling(exchanged_boxes_size(j_box)%line(i_box)%&
                    exchanged%get_frequency_ratio() * real(nums_candidates(i_couple), DP))
                nums_candidates(i_couple) = merge(nums_candidates(i_couple), 1, &
                    nums_candidates(i_couple) > 0)
            else
                nums_candidates(i_couple) = 0
            end if
        end do

        call selector%reset(nums_candidates)
    end subroutine reset_boxes_volume_exchange

end module procedures_selectors_resetters
