module procedures_energies_writers_factory

use types_string_wrapper, only: String_Wrapper
use classes_external_field, only: Abstract_External_Field
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_total_moments_factory, only: set_are_dipolar
use classes_pair_potential, only: Pair_Potential_Wrapper, Pair_Potentials_Line
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

interface create
    module procedure :: create_line
    module procedure :: create_element
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_line
end interface destroy

contains

    subroutine create_line(energies, paths, external_fields, wall_pairs, components, short_pairs, &
        visit_energies)
        type(Concrete_Energies_Writers), allocatable, intent(out) :: energies(:)
        type(String_Wrapper), intent(in) :: paths(:)
        class(Abstract_External_Field), intent(in) :: external_fields(:)
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        type(Pair_Potentials_Line), intent(in) :: short_pairs(:)
        logical, intent(in) :: visit_energies

        logical :: are_dipolar(size(components, 1), size(components, 2))
        integer :: i_box

        call set_are_dipolar(are_dipolar, components)
        allocate(energies(size(paths)))
        do i_box = 1, size(energies)
            call create(energies(i_box), paths(i_box)%string, external_fields(i_box), wall_pairs, &
                components(:, i_box), are_dipolar(:, i_box), short_pairs, visit_energies)
        end do
    end subroutine create_line

    subroutine destroy_line(energies)
        type(Concrete_Energies_Writers), allocatable, intent(inout) :: energies(:)

        integer :: i_box

        if (allocated(energies)) then
            do i_box = size(energies), 1, -1
                call destroy(energies(i_box))
            end do
            deallocate(energies)
        end if
    end subroutine destroy_line

    subroutine create_element(energies, path, external_field, wall_pairs, components, are_dipolar, &
        short_pairs, visit_energies)
        type(Concrete_Energies_Writers), intent(out) :: energies
        character(len=*), intent(in) :: path
        class(Abstract_External_Field), intent(in) :: external_field
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)
        type(Pair_Potentials_Line), intent(in) :: short_pairs(:)
        logical, intent(in) :: visit_energies

        call line_writer_create(energies%walls_energies, path//"walls_energies.out", wall_pairs, &
            visit_energies)
        call line_writer_create(energies%field_energies, path//"field_energies.out", &
            external_field, components, visit_energies)
        call triangle_writer_create(energies%short_energies, path//"short_energies.out", &
            short_pairs, visit_energies)
        call triangle_writer_create(energies%dipolar_energies, path//"dipolar_energies.out", &
            are_dipolar, visit_energies)
        call real_writer_create(energies%dipolar_shared_energy, path//"dipolar_shared_energy.out", &
            any(are_dipolar) .and. visit_energies)
    end subroutine create_element

    subroutine destroy_element(energies)
        type(Concrete_Energies_Writers), intent(inout) :: energies

        call real_writer_destroy(energies%dipolar_shared_energy)
        call triangle_writer_destroy(energies%dipolar_energies)
        call triangle_writer_destroy(energies%short_energies)
        call line_writer_destroy(energies%field_energies)
        call line_writer_destroy(energies%walls_energies)
    end subroutine destroy_element

end module procedures_energies_writers_factory
