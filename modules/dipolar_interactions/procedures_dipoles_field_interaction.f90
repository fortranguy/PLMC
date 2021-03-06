module procedures_dipoles_field_interaction

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_external_field, only: Abstract_External_Field
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments
use types_particle_wrapper, only: Concrete_Particle

implicit none

private
public :: visit_component, visit_translation, visit_rotation, visit_add, visit_remove

contains

    pure real(DP) function visit_component(external_field, positions, dipole_moments) &
        result(energy)
        class(Abstract_External_Field), intent(in) :: external_field
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Component_Dipole_Moments), intent(in) :: dipole_moments

        integer :: i_particle

        energy = 0
        do i_particle = 1, dipole_moments%get_num()
            energy = energy - dot_product(dipole_moments%get(i_particle), external_field%&
                get(positions%get(i_particle)))
        end do
    end function visit_component

    pure real(DP) function visit_translation(external_field, new_position, old) result(delta_energy)
        class(Abstract_External_Field), intent(in) :: external_field
        real(DP), dimension(:), intent(in) :: new_position
        type(Concrete_Particle), intent(in) :: old

        delta_energy = -dot_product(old%dipole_moment, external_field%get(new_position) - &
            external_field%get(old%position))
    end function visit_translation

    pure real(DP) function visit_rotation(external_field, new_dipole_moment, old) &
        result(delta_energy)
        class(Abstract_External_Field), intent(in) :: external_field
        real(DP), dimension(:), intent(in) :: new_dipole_moment
        type(Concrete_Particle), intent(in) :: old

        delta_energy = -dot_product(new_dipole_moment - old%dipole_moment, external_field%&
            get(old%position))
    end function visit_rotation

    pure real(DP) function visit_add(external_field, particle) result(delta_energy)
        class(Abstract_External_Field), intent(in) :: external_field
        type(Concrete_Particle), intent(in) :: particle

        delta_energy = visit_exchange(external_field, particle, +1._DP)
    end function visit_add

    pure real(DP) function visit_remove(external_field, particle) result(delta_energy)
        class(Abstract_External_Field), intent(in) :: external_field
        type(Concrete_Particle), intent(in) :: particle

        delta_energy = visit_exchange(external_field, particle, -1._DP)
    end function visit_remove

    pure real(DP) function visit_exchange(external_field, particle, signed) result(delta_energy)
        class(Abstract_External_Field), intent(in) :: external_field
        type(Concrete_Particle), intent(in) :: particle
        real(DP), intent(in) :: signed

        delta_energy = -signed * dot_product(particle%dipole_moment, external_field%get(particle%&
            position))
    end function visit_exchange

end module procedures_dipoles_field_interaction
