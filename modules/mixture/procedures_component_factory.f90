module procedures_component_factory

use json_module, only: json_file
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_coordinates_factory, only: coordinates_create => create, &
    coordinates_destroy => destroy
use procedures_composition_factory, only: composition_create => create, &
    composition_destroy => destroy
use types_component_wrapper, only: Component_Wrapper
use procedures_property_inquirers, only: component_is_dipolar, component_can_exchange

implicit none

private
public :: create, destroy

contains

    subroutine create(component, exists, periodic_box, accessible_domain, generating_data, prefix)
        type(Component_Wrapper), intent(out) :: component
        logical, intent(in) :: exists
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Parallelepiped_Domain), intent(in) :: accessible_domain
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        logical :: is_dipolar, can_exchange

        call composition_create(component%number, exists)
        call coordinates_create(component%positions, exists, periodic_box, component%number)
        is_dipolar = exists .and. component_is_dipolar(generating_data, prefix)
        call coordinates_create(component%orientations, component%number, is_dipolar)
        call coordinates_create(component%dipole_moments, component%orientations, is_dipolar, &
            generating_data, prefix)
        can_exchange = exists .and. component_can_exchange(generating_data, prefix)
        call composition_create(component%chemical_potential, can_exchange, generating_data, prefix)
        call composition_create(component%average_number, accessible_domain, component%number, &
            component%chemical_potential)
    end subroutine create

    subroutine destroy(component)
        type(Component_Wrapper), intent(inout) :: component

        call composition_destroy(component%average_number)
        call composition_destroy(component%chemical_potential)
        call coordinates_destroy(component%dipole_moments)
        call coordinates_destroy(component%orientations)
        call coordinates_destroy(component%positions)
        call composition_destroy(component%number)
    end subroutine destroy

end module procedures_component_factory
