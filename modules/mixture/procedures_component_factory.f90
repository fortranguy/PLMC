module procedures_component_factory

use json_module, only: json_file
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_composition_factory, only: composition_create => create, &
    composition_destroy => destroy
use procedures_coordinates_factory, only: coordinates_create => create, &
    coordinates_destroy => destroy
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only: component_is_dipolar

implicit none

private
public :: create, destroy

contains

    subroutine create(component, periodic_box, exists, can_exchange, generating_data, prefix)
        type(Component_Wrapper), intent(out) :: component
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: exists, can_exchange
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        logical :: is_dipolar

        call composition_create(component%num_particles, exists)
        call coordinates_create(component%positions, periodic_box, component%num_particles, exists)
        is_dipolar = exists .and. component_is_dipolar(generating_data, prefix)
        call coordinates_create(component%orientations, component%num_particles, is_dipolar)
        call coordinates_create(component%dipole_moments, component%orientations, is_dipolar, &
            generating_data, prefix)
        call composition_create(component%chemical_potential, can_exchange, generating_data, prefix)
    end subroutine create

    subroutine destroy(component)
        type(Component_Wrapper), intent(inout) :: component

        call composition_destroy(component%chemical_potential)
        call coordinates_destroy(component%dipole_moments)
        call coordinates_destroy(component%orientations)
        call coordinates_destroy(component%positions)
        call composition_destroy(component%num_particles)
    end subroutine destroy

end module procedures_component_factory
