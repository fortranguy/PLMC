module procedures_dipoles_field_interaction

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_external_field, only: Abstract_External_Field
use types_temporary_particle, only: Concrete_Temporary_Particle

implicit none

private

contains

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
