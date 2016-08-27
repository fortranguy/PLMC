module procedures_energies_writers_factory

use classes_external_field, only: Abstract_External_Field
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_total_moment_factory, only: set_are_dipolar
use types_pair_potential_wrapper, only: Pair_Potential_Wrapper, Pair_Potentials_Line
use procedures_real_writer_factory, only: real_writer_create => create, &
    real_writer_destroy => destroy
use procedures_line_writer_factory, only: line_writer_create => create, &
    line_writer_destroy => destroy
use procedures_triangle_writer_factory, only: triangle_writer_create => create, &
    triangle_writer_destroy => destroy
use types_energies_writers, only: Concrete_Energies_Writers

implicit none

private
public :: create, destroy

contains

    subroutine create(energies, external_field, wall_pairs, components, short_pairs, visit_energies)
        type(Concrete_Energies_Writers), intent(out) :: energies
        class(Abstract_External_Field), intent(in) :: external_field
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Pair_Potentials_Line), intent(in) :: short_pairs(:)
        logical, intent(in) :: visit_energies

        logical, dimension(size(components)) :: are_dipolar

        call line_writer_create(energies%walls_energies, wall_pairs, visit_energies, &
            "walls_energies.out")
        call line_writer_create(energies%field_energies, external_field, components, &
            visit_energies, "field_energies.out")
        call triangle_writer_create(energies%short_energies, short_pairs, visit_energies, &
            "short_energies.out")
        call set_are_dipolar(are_dipolar, components)
        call triangle_writer_create(energies%dipolar_energies, are_dipolar, visit_energies, &
            "dipolar_energies.out")
        call real_writer_create(energies%dipolar_mixture_energy, any(are_dipolar) .and. &
            visit_energies, "dipolar_mixture_energy.out")
    end subroutine create

    subroutine destroy(energies)
        type(Concrete_Energies_Writers), intent(inout) :: energies

        call real_writer_destroy(energies%dipolar_mixture_energy)
        call triangle_writer_destroy(energies%dipolar_energies)
        call triangle_writer_destroy(energies%short_energies)
        call line_writer_destroy(energies%field_energies)
        call line_writer_destroy(energies%walls_energies)
    end subroutine destroy

end module procedures_energies_writers_factory
