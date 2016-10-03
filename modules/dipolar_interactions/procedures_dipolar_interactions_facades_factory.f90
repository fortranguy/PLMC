module procedures_dipolar_interactions_facades_factory

use json_module, only: json_file
use procedures_errors, only: error_exit
use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_inquirers, only: periodicity_is_xyz, box_size_can_change
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_total_moments_factory, only: set_are_dipolar
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use classes_dipolar_interactions_facade, only: Abstract_Dipolar_Interactions_Facade, &
    Scalable_Dipolar_Interactions_Facade, Unscalable_Dipolar_Interactions_Facade, &
    Null_Dipolar_Interactions_Facade
use procedures_exploration_inquirers, only: property_measure_pressure => measure_pressure

implicit none

private
public :: create, destroy

contains

    subroutine create(facades, environment, components, dipolar_interactions_dynamic, &
        dipolar_interactions_static, exploring_data, volume_change_prefix)
        class(Abstract_Dipolar_Interactions_Facade), allocatable, intent(out) :: facades(:)
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:, :)
        type(Dipolar_Interactions_Dynamic_Wrapper), intent(in) :: dipolar_interactions_dynamic
        type(Dipolar_Interactions_Static_Wrapper), intent(in) :: dipolar_interactions_static(:)
        type(json_file), optional, intent(inout) :: exploring_data
        character(len=*), optional, intent(in) :: volume_change_prefix

        integer :: i_box
        logical :: are_dipolar(size(components, 1), size(components, 2)), measure_pressure

        if (present(exploring_data) .and. present(volume_change_prefix)) then
            measure_pressure = property_measure_pressure(exploring_data, volume_change_prefix)
        else
            measure_pressure = .false.
        end if

        do i_box = 1, size(are_dipolar, 2)
            call set_are_dipolar(are_dipolar(:, i_box), components(:, i_box))
        end do

        if ((box_size_can_change(environment%beta_pressure) .or. measure_pressure) .and. &
            any(are_dipolar)) then
            if (all(periodicity_is_xyz(environment%periodic_boxes))) then
                allocate(Scalable_Dipolar_Interactions_Facade :: &
                    facades(size(dipolar_interactions_static)))
            else
                allocate(Unscalable_Dipolar_Interactions_Facade :: &
                    facades(size(dipolar_interactions_static)))
            end if
        else
            allocate(Null_Dipolar_Interactions_Facade :: facades(size(dipolar_interactions_static)))
        end if

        select type (facades)
            type is (Scalable_Dipolar_Interactions_Facade)
                do i_box = 1, size(facades)
                    call facades(i_box)%construct(components(:, i_box), &
                        dipolar_interactions_dynamic, dipolar_interactions_static(i_box))
                end do
            type is (Unscalable_Dipolar_Interactions_Facade)
                do i_box = 1, size(facades)
                    call facades(i_box)%construct(environment%periodic_boxes(i_box), &
                        components(:, i_box), dipolar_interactions_dynamic, &
                        dipolar_interactions_static(i_box))
                end do
            type is (Null_Dipolar_Interactions_Facade)
            class default
                call error_exit("procedures_dipolar_interactions_facades_factory: create: "//&
                    "facades type is unknown.")
        end select
    end subroutine create

    subroutine destroy(facades)
        class(Abstract_Dipolar_Interactions_Facade), allocatable, intent(inout) :: facades(:)

        integer :: i_box

        if (allocated(facades)) then
            do i_box = size(facades), 1, -1
                call facades(i_box)%destroy()
            end do
            deallocate(facades)
        end if
    end subroutine destroy

end module procedures_dipolar_interactions_facades_factory
