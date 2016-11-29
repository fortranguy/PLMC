module classes_structure_factor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_particle_wrapper, only: Concrete_Particle

implicit none

private

    type, abstract, public :: Abstract_Structure_Factor
    contains
        procedure(Abstract_is_dipolar), deferred :: is_dipolar
        procedure(Abstract_update_translation), deferred :: update_translation
        procedure(Abstract_update_transmutation), deferred :: update_transmutation
        procedure(Abstract_update_rotation), deferred :: update_rotation
        procedure(Abstract_update_add), deferred :: update_add
        procedure(Abstract_update_remove), deferred :: update_remove
        procedure(Abstract_update_switch), deferred :: update_switch
    end type Abstract_Structure_Factor

    abstract interface

        pure logical function Abstract_is_dipolar(this, i_component) result(is_dipolar)
        import :: Abstract_Structure_Factor
            class(Abstract_Structure_Factor), intent(in) :: this
            integer, intent(in) :: i_component
        end function Abstract_is_dipolar

        pure subroutine Abstract_update_translation(this, i_component, new_position, old)
        import :: DP, Concrete_Particle, Abstract_Structure_Factor
            class(Abstract_Structure_Factor), intent(inout) :: this
            integer, intent(in) :: i_component
            real(DP), intent(in) :: new_position(:)
            type(Concrete_Particle), intent(in) :: old
        end subroutine Abstract_update_translation

        pure subroutine Abstract_update_transmutation(this, ij_components, new_dipole_moment, old)
        import :: DP, Concrete_Particle, Abstract_Structure_Factor
            class(Abstract_Structure_Factor), intent(inout) :: this
            integer, intent(in) :: ij_components(:)
            real(DP), intent(in) :: new_dipole_moment(:)
            type(Concrete_Particle), intent(in) :: old
        end subroutine Abstract_update_transmutation

        pure subroutine Abstract_update_rotation(this, i_component, new_dipole_moment, old)
        import :: DP, Concrete_Particle, Abstract_Structure_Factor
            class(Abstract_Structure_Factor), intent(inout) :: this
            integer, intent(in) :: i_component
            real(DP), intent(in) :: new_dipole_moment(:)
            type(Concrete_Particle), intent(in) :: old
        end subroutine Abstract_update_rotation

        pure subroutine Abstract_update_add(this, i_component, particle)
        import :: Concrete_Particle, Abstract_Structure_Factor
            class(Abstract_Structure_Factor), intent(inout) :: this
            integer, intent(in) :: i_component
            type(Concrete_Particle), intent(in) :: particle
        end subroutine Abstract_update_add

        pure subroutine Abstract_update_remove(this, i_component, particle)
        import :: Concrete_Particle, Abstract_Structure_Factor
            class(Abstract_Structure_Factor), intent(inout) :: this
            integer, intent(in) :: i_component
            type(Concrete_Particle), intent(in) :: particle
        end subroutine Abstract_update_remove

        pure subroutine Abstract_update_switch(this, ij_components, particles)
        import :: Concrete_Particle, Abstract_Structure_Factor
            class(Abstract_Structure_Factor), intent(inout) :: this
            integer, intent(in) :: ij_components(:)
            type(Concrete_Particle), intent(in) :: particles(:)
        end subroutine Abstract_update_switch

    end interface

end module classes_structure_factor
