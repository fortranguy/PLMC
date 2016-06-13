module classes_structure_visitor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_temporary_particle, only: Concrete_Temporary_Particle

implicit none

private

    type, abstract, public :: Abstract_Structure_Visitor
    contains
        procedure(Abstract_visit), deferred :: visit
        procedure(Abstract_visit_translation), deferred :: visit_translation
        procedure(Abstract_visit_transmutation), deferred :: visit_transmutation
        procedure(Abstract_visit_rotation), deferred :: visit_rotation
        procedure(Abstract_visit_add), deferred :: visit_add
        procedure(Abstract_visit_remove), deferred :: visit_remove
        procedure(Abstract_visit_switch), deferred :: visit_switch
    end type Abstract_Structure_Visitor

    abstract interface

        pure real(DP) function Abstract_visit(this)
        import :: DP, Abstract_Structure_Visitor
            class(Abstract_Structure_Visitor), intent(in) :: this
        end function Abstract_visit

        pure real(DP) function Abstract_visit_translation(this, i_component, new_position, old)
        import :: DP, Concrete_Temporary_Particle, Abstract_Structure_Visitor
            class(Abstract_Structure_Visitor), intent(in) :: this
            integer, intent(in) :: i_component
            real(DP), intent(in) :: new_position(:)
            type(Concrete_Temporary_Particle), intent(in) :: old
        end function Abstract_visit_translation

        pure real(DP) function Abstract_visit_transmutation(this, ij_components,new_dipolar_moment,&
            old)
        import :: DP, Concrete_Temporary_Particle, Abstract_Structure_Visitor
            class(Abstract_Structure_Visitor), intent(in) :: this
            integer, intent(in) :: ij_components(:)
            real(DP), intent(in) :: new_dipolar_moment(:)
            type(Concrete_Temporary_Particle), intent(in) :: old
        end function Abstract_visit_transmutation

        pure real(DP) function Abstract_visit_rotation(this, i_component, new_dipolar_moment, old)
        import :: DP, Concrete_Temporary_Particle, Abstract_Structure_Visitor
            class(Abstract_Structure_Visitor), intent(in) :: this
            integer, intent(in) :: i_component
            real(DP), intent(in) :: new_dipolar_moment(:)
            type(Concrete_Temporary_Particle), intent(in) :: old
        end function Abstract_visit_rotation

        pure real(DP) function Abstract_visit_add(this, i_component, particle)
        import :: DP, Concrete_Temporary_Particle, Abstract_Structure_Visitor
        class(Abstract_Structure_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle
        end function Abstract_visit_add

        pure real(DP) function Abstract_visit_remove(this, i_component, particle)
        import :: DP, Concrete_Temporary_Particle, Abstract_Structure_Visitor
        class(Abstract_Structure_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle
        end function Abstract_visit_remove

        pure real(DP) function Abstract_visit_switch(this, ij_components, particles)
        import :: DP, Concrete_Temporary_Particle, Abstract_Structure_Visitor
            class(Abstract_Structure_Visitor), intent(in) :: this
            integer, intent(in) :: ij_components(:)
            type(Concrete_Temporary_Particle), intent(in) :: particles(:)
        end function Abstract_visit_switch

    end interface

end module classes_structure_visitor
