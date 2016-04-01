module procedures_dipoles_field_interaction

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_external_field, only: Abstract_External_Field
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_temporary_particle, only: Concrete_Temporary_Particle

implicit none

private
public :: dipoles_field_visit, dipoles_field_visit_move, dipoles_field_visit_rotation

contains

    pure real(DP) function dipoles_field_visit(external_field, component_dipolar_moments, &
        component_positions) result(energy)
        class(Abstract_External_Field), intent(in) :: external_field
        class(Abstract_Component_Dipolar_Moments), intent(in) :: component_dipolar_moments
        class(Abstract_Component_Coordinates), intent(in) :: component_positions

        integer :: i_particle

        energy = 0
        do i_particle = 1, component_dipolar_moments%get_num()
            energy = energy -dot_product(component_dipolar_moments%get(i_particle), external_field%&
                get(component_positions%get(i_particle)))
        end do
    end function dipoles_field_visit

    pure real(DP) function dipoles_field_visit_move(external_field, new_position, old) &
        result(delta_energy)
        class(Abstract_External_Field), intent(in) :: external_field
        real(DP), dimension(:), intent(in) :: new_position
        type(Concrete_Temporary_Particle), intent(in) :: old

        delta_energy = -dot_product(old%dipolar_moment, external_field%get(new_position) - &
            external_field%get(old%position))
    end function dipoles_field_visit_move

    pure real(DP) function dipoles_field_visit_rotation(external_field, new_dipolar_moment, old) &
        result(delta_energy)
        class(Abstract_External_Field), intent(in) :: external_field
        real(DP), dimension(:), intent(in) :: new_dipolar_moment
        type(Concrete_Temporary_Particle), intent(in) :: old

        delta_energy = -dot_product(new_dipolar_moment - old%dipolar_moment, external_field%&
            get(old%position))
    end function dipoles_field_visit_rotation

end module procedures_dipoles_field_interaction
