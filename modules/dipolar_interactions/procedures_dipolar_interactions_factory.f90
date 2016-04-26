module procedures_dipolar_interactions_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: warning_continue
use classes_permittivity, only: Abstract_Permittivity
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use procedures_mixture_total_moment_factory, only: set_are_dipolar
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use procedures_des_convergence_parameter_factory, only: des_convergence_parameter_create => create,&
    des_convergence_parameter_destroy => destroy
use procedures_des_real_factory, only: des_real_create => create, des_real_destroy => destroy
use procedures_des_reci_factory, only: des_reci_create => create, des_reci_destroy => destroy
use procedures_des_self_factory, only: des_self_create => create, des_self_destroy => destroy
use procedures_des_surf_mixture_factory, only: des_surf_mixture_create => create, &
    des_surf_mixture_destroy => destroy
use procedures_dlc_factory, only: dlc_create => create, dlc_destroy => destroy
use procedures_property_inquirers, only: use_permittivity, use_reciprocal_lattice

implicit none

private
public :: dipolar_interactions_create, dipolar_interactions_destroy

contains

    subroutine dipolar_interactions_create(dipolar_interactions, environment, mixture, input_data, &
        prefix)
        type(Dipolar_Interactions_Wrapper), intent(out) :: dipolar_interactions
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        logical :: are_dipolar(size(mixture%components))

        call set_are_dipolar(are_dipolar, mixture%components)
        call check_consistency(environment%reciprocal_lattice, environment%permittivity, &
            any(are_dipolar))
        call des_convergence_parameter_create(dipolar_interactions%alpha, environment%periodic_box,&
            any(are_dipolar), input_data, prefix)

        call des_real_create(dipolar_interactions%real_visitor, environment%periodic_box, &
            any(are_dipolar))
        call des_real_create(dipolar_interactions%real_pair, environment%permittivity, mixture%&
            components_min_distances, any(are_dipolar), dipolar_interactions%alpha, input_data, &
            prefix//"Real.")
        call des_real_create(dipolar_interactions%real_components, environment%periodic_box, &
            mixture%components, are_dipolar, dipolar_interactions%real_pair)

        call des_reci_create(dipolar_interactions%reci_weight, environment, any(are_dipolar), &
            dipolar_interactions%alpha)
        call des_reci_create(dipolar_interactions%reci_structure, environment, mixture%components, &
            are_dipolar)
        call des_reci_create(dipolar_interactions%reci_visitor, environment, &
            dipolar_interactions%reci_weight, dipolar_interactions%reci_structure)

        call des_self_create(dipolar_interactions%self_components, environment%permittivity, &
            mixture%components, are_dipolar, dipolar_interactions%alpha)

        call des_surf_mixture_create(dipolar_interactions%surf_mixture, environment%periodic_box, &
            environment%permittivity, mixture%total_moment)

        call dlc_create(dipolar_interactions%dlc_weight, environment, any(are_dipolar))
        call dlc_create(dipolar_interactions%dlc_structures, environment, mixture%components, &
            are_dipolar)
        call dlc_create(dipolar_interactions%dlc_visitor, environment, dipolar_interactions%&
            dlc_weight, dipolar_interactions%dlc_structures)
    end subroutine dipolar_interactions_create

    subroutine dipolar_interactions_destroy(dipolar_interactions)
        type(Dipolar_Interactions_Wrapper), intent(inout) :: dipolar_interactions

        call dlc_destroy(dipolar_interactions%dlc_visitor)
        call dlc_destroy(dipolar_interactions%dlc_structures)
        call dlc_destroy(dipolar_interactions%dlc_weight)

        call des_surf_mixture_destroy(dipolar_interactions%surf_mixture)

        call des_self_destroy(dipolar_interactions%self_components)

        call des_reci_destroy(dipolar_interactions%reci_visitor)
        call des_reci_destroy(dipolar_interactions%reci_structure)
        call des_reci_destroy(dipolar_interactions%reci_weight)

        call des_real_destroy(dipolar_interactions%real_components)
        call des_real_destroy(dipolar_interactions%real_pair)
        call des_real_destroy(dipolar_interactions%real_visitor)

        call des_convergence_parameter_destroy(dipolar_interactions%alpha)
    end subroutine dipolar_interactions_destroy

    subroutine check_consistency(reciprocal_lattice, permittivity, dipoles_exist)
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Permittivity), intent(in) :: permittivity
        logical, intent(in) :: dipoles_exist

        if (dipoles_exist) then
            if (.not.use_reciprocal_lattice(reciprocal_lattice)) then
                call warning_continue("dipolar_interactions_check: "//&
                    "dipoles exist but reciprocal_lattice unused.")
            end if
            if (.not.use_permittivity(permittivity)) then
                call warning_continue("dipolar_interactions_check: "//&
                    "dipoles exist but permittivity unused.")
            end if
        end if
    end subroutine check_consistency

end module procedures_dipolar_interactions_factory
