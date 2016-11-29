module procedures_mixture_properties

use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only: component_exists, component_has_positions, &
    component_has_orientations, component_can_exchange

implicit none

private
public :: set_exist, set_nums_particles, set_have_positions, set_have_orientations, set_can_exchange

contains

    subroutine set_exist(components_exist, components)
        logical, intent(inout) :: components_exist(:, :)
        type(Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box, i_component

        do i_box = 1, size(components_exist, 2)
            do i_component = 1, size(components_exist, 1)
                components_exist(i_component, i_box) = &
                    component_exists(components(i_component, i_box)%num_particles)
            end do
        end do
    end subroutine set_exist

    subroutine set_nums_particles(nums_particles, components)
        integer, intent(inout) :: nums_particles(:, :)
        type(Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box, i_component

        do i_box = 1, size(nums_particles, 2)
            do i_component = 1, size(nums_particles, 1)
                nums_particles(i_component, i_box) = components(i_component, i_box)%num_particles%&
                    get()
            end do
        end do
    end subroutine set_nums_particles

    subroutine set_have_positions(have_positions, components)
        logical, intent(inout) :: have_positions(:, :)
        type(Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box, i_component

        do i_box = 1, size(have_positions, 2)
            do i_component = 1, size(have_positions, 1)
                have_positions(i_component, i_box) = &
                    component_has_positions(components(i_component, i_box)%positions)
            end do
        end do
    end subroutine set_have_positions

    subroutine set_have_orientations(have_orientations, components)
        logical, intent(inout) :: have_orientations(:, :)
        type(Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box, i_component

        do i_box = 1, size(have_orientations, 2)
            do i_component = 1, size(have_orientations, 1)
                have_orientations(i_component, i_box) = &
                    component_has_orientations(components(i_component, i_box)%orientations)
            end do
        end do
    end subroutine set_have_orientations

    subroutine set_can_exchange(can_exchange, components)
        logical, intent(inout) :: can_exchange(:, :)
        type(Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box, i_component

        do i_box = 1, size(can_exchange, 2)
            do i_component = 1, size(can_exchange, 1)
                can_exchange(i_component, i_box) = &
                    component_can_exchange(components(i_component, i_box)%chemical_potential)
            end do
        end do
    end subroutine set_can_exchange

end module procedures_mixture_properties
