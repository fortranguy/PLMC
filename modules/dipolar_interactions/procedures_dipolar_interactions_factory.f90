module procedures_dipolar_interactions_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: warning_continue
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_permittivity, only: Abstract_Permittivity
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use procedures_mixture_factory, only: mixture_set
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter, &
    Concrete_DES_Convergence_Parameter, Null_DES_Convergence_Parameter
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use procedures_property_inquirers, only: use_permittivity, use_reciprocal_lattice
use procedures_des_real_factory, only: des_real_create, des_real_destroy
use procedures_des_reci_factory, only: des_reci_create, des_reci_destroy
use procedures_des_self_factory, only: des_self_create, des_self_destroy
use procedures_des_surf_factory, only: des_surf_create, des_surf_destroy
use procedures_dlc_factory, only: dlc_create, dlc_destroy

implicit none

private
public :: dipolar_interactions_create, dipolar_interactions_destroy

interface dipolar_interactions_create
    module procedure :: create_all
    module procedure :: create_alpha
end interface dipolar_interactions_create

interface dipolar_interactions_destroy
    module procedure :: destroy_alpha
    module procedure :: destroy_all
end interface dipolar_interactions_destroy

interface dipolar_interactions_check
    module procedure :: check_consistency
end interface dipolar_interactions_check

contains

    subroutine create_all(dipolar_interactions, environment, mixture, input_data, prefix)
        type(Dipolar_Interactions_Wrapper), intent(out) :: dipolar_interactions
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        logical :: are_dipolar(size(mixture%components))

        call mixture_set(are_dipolar, mixture%components)
        call dipolar_interactions_check(environment%reciprocal_lattice, environment%permittivity, &
            any(are_dipolar))
        call dipolar_interactions_create(dipolar_interactions%alpha, environment%periodic_box, &
            any(are_dipolar), input_data, prefix)

        call des_real_create(dipolar_interactions%real_visitor, environment%periodic_box, &
            any(are_dipolar))
        call des_real_create(dipolar_interactions%real_pairs, environment%permittivity, mixture%&
            components_min_distances, are_dipolar, dipolar_interactions%alpha, input_data, &
            prefix//"Real.")
        call des_real_create(dipolar_interactions%real_components, environment%periodic_box, &
            mixture%components, dipolar_interactions%real_pairs)

        call des_reci_create(dipolar_interactions%reci_weight, environment, any(are_dipolar), &
            dipolar_interactions%alpha)
        call des_reci_create(dipolar_interactions%reci_structure, environment, mixture%components, &
            are_dipolar)
        call des_reci_create(dipolar_interactions%reci_visitor, environment, &
            dipolar_interactions%reci_weight, dipolar_interactions%reci_structure)

        call des_self_create(dipolar_interactions%self_components, environment%permittivity, &
            mixture%components, are_dipolar, dipolar_interactions%alpha)

        call des_surf_create(dipolar_interactions%surf_mixture, environment%periodic_box, &
            environment%permittivity, mixture%total_moment)

        call dlc_create(dipolar_interactions%dlc_weight, environment, any(are_dipolar))
        call dlc_create(dipolar_interactions%dlc_structures, environment, mixture%components, &
            are_dipolar)
        call dlc_create(dipolar_interactions%dlc_visitor, environment, dipolar_interactions%&
            dlc_weight, dipolar_interactions%dlc_structures)
    end subroutine create_all

    subroutine destroy_all(dipolar_interactions)
        type(Dipolar_Interactions_Wrapper), intent(inout) :: dipolar_interactions

        call dlc_destroy(dipolar_interactions%dlc_visitor)
        call dlc_destroy(dipolar_interactions%dlc_structures)
        call dlc_destroy(dipolar_interactions%dlc_weight)

        call des_surf_destroy(dipolar_interactions%surf_mixture)

        call des_self_destroy(dipolar_interactions%self_components)

        call des_reci_destroy(dipolar_interactions%reci_visitor)
        call des_reci_destroy(dipolar_interactions%reci_structure)
        call des_reci_destroy(dipolar_interactions%reci_weight)

        call des_real_destroy(dipolar_interactions%real_components)
        call des_real_destroy(dipolar_interactions%real_pairs)
        call des_real_destroy(dipolar_interactions%real_visitor)

        call dipolar_interactions_destroy(dipolar_interactions%alpha)
    end subroutine destroy_all

    subroutine create_alpha(alpha, periodic_box, dipoles_exist, input_data, prefix)
        class(Abstract_DES_Convergence_Parameter), allocatable, intent(out) :: alpha
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: dipoles_exist
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: alpha_x_box

        if (dipoles_exist) then
            data_field = prefix//"alpha times box edge"
            call input_data%get(data_field, alpha_x_box, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_DES_Convergence_Parameter :: alpha)
        else
            allocate(Null_DES_Convergence_Parameter :: alpha)
        end if
        call alpha%construct(periodic_box, alpha_x_box)
    end subroutine create_alpha

    subroutine destroy_alpha(alpha)
        class(Abstract_DES_Convergence_Parameter), allocatable, intent(inout) :: alpha

        if (allocated(alpha)) then
            call alpha%destroy()
            deallocate(alpha)
        end if
    end subroutine destroy_alpha

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
