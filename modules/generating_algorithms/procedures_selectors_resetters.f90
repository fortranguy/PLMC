module procedures_selectors_resetters

use classes_tower_sampler, only: Abstract_Tower_Sampler
use types_component_wrapper, only: Component_Wrapper

implicit none

private
public :: reset

interface reset
    module procedure :: reset_one_particle
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

end module procedures_selectors_resetters
