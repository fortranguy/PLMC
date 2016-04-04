module procedures_dipoles_field_interaction

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_external_field, only: Abstract_External_Field
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_temporary_particle, only: Concrete_Temporary_Particle

implicit none

private
public :: dipoles_field_visit_component, dipoles_field_visit_move, dipoles_field_visit_rotation

contains

    pure real(DP) function dipoles_field_visit_component(external_field, positions, &
        dipolar_moments) result(energy)
        class(Abstract_External_Field), intent(in) :: external_field
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Component_Dipolar_Moments), intent(in) :: dipolar_moments

        integer :: i_particle

        energy = 0
        do i_particle = 1, dipolar_moments%get_num()
            energy = energy - dot_product(dipolar_moments%get(i_particle), external_field%&
                get(positions%get(i_particle)))
        end do
    end function dipoles_field_visit_component

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

    pure real(DP) function dipoles_field_visit_switch(external_field, particles) &
        result(delta_energy)
        class(Abstract_External_Field), intent(in) :: external_field
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        delta_energy = dot_product(particles(2)%dipolar_moment - particles(1)%dipolar_moment, &
            external_field%get(particles(2)%position) - external_field%get(particles(1)%position))
    end function dipoles_field_visit_switch

end module procedures_dipoles_field_interaction