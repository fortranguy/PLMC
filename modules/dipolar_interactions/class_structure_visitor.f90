module class_structure_visitor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_temporary_particle, only: Concrete_Temporary_Particle

implicit none

private

    type, abstract, public :: Abstract_Structure_Visitor
    contains
        procedure(Abstract_visit), deferred :: visit
        procedure(Abstract_visit_move), deferred :: visit_move
        procedure(Abstract_visit_rotation), deferred :: visit_rotation
        procedure(Abstract_visit_switch), deferred :: visit_switch
    end type Abstract_Structure_Visitor

    abstract interface

        pure real(DP) function Abstract_visit(this) result(energy)
        import :: DP, Abstract_Structure_Visitor
            class(Abstract_Structure_Visitor), intent(in) :: this
        end function Abstract_visit

        pure real(DP) function Abstract_visit_move(this, i_component, new_position, old) &
            result(delta_energy)
        import :: DP, Concrete_Temporary_Particle, Abstract_Structure_Visitor
            class(Abstract_Structure_Visitor), intent(in) :: this
            integer, intent(in) :: i_component
            real(DP), intent(in) :: new_position(:)
            type(Concrete_Temporary_Particle), intent(in) :: old
        end function Abstract_visit_move

        pure real(DP) function Abstract_visit_rotation(this, i_component, new_dipolar_moment, old) &
            result(delta_energy)
        import :: DP, Concrete_Temporary_Particle, Abstract_Structure_Visitor
            class(Abstract_Structure_Visitor), intent(in) :: this
            integer, intent(in) :: i_component
            real(DP), intent(in) :: new_dipolar_moment(:)
            type(Concrete_Temporary_Particle), intent(in) :: old
        end function Abstract_visit_rotation

        pure real(DP) function Abstract_visit_switch(this, ij_components, particles) &
            result(delta_energy)
        import :: DP, Concrete_Temporary_Particle, Abstract_Structure_Visitor
            class(Abstract_Structure_Visitor), intent(in) :: this
            integer, intent(in) :: ij_components(:)
            type(Concrete_Temporary_Particle), intent(in) :: particles(:)
        end function Abstract_visit_switch

    end interface

end module class_structure_visitor
