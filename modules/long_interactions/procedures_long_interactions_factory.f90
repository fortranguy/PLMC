module procedures_long_interactions_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: warning_continue
use procedures_checks, only: check_data_found
use class_periodic_box, only: Abstract_Periodic_Box
use class_permittivity, only: Abstract_Permittivity
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use class_ewald_convergence_parameter, only: Abstract_Ewald_Convergence_Parameter, &
    Concrete_Ewald_Convergence_Parameter, Null_Ewald_Convergence_Parameter
use types_long_interactions_wrapper, only: Long_Interactions_Wrapper
use procedures_property_inquirers, only: use_permittivity, use_reciprocal_lattice, &
    component_is_dipolar
use procedures_ewald_real_factory, only: ewald_real_create, ewald_real_destroy
use procedures_ewald_reci_factory, only: ewald_reci_create, ewald_reci_destroy

implicit none

private
public :: long_interactions_create, long_interactions_set, long_interactions_destroy

interface long_interactions_create
    module procedure :: create_all
    module procedure :: create_alpha
end interface long_interactions_create

interface long_interactions_set
    module procedure :: set_are_dipolar
end interface long_interactions_set

interface long_interactions_check
    module procedure :: check_consistency
end interface long_interactions_check

interface long_interactions_destroy
    module procedure :: destroy_alpha
    module procedure :: destroy_all
end interface long_interactions_destroy

contains

    subroutine create_all(long_interactions, environment, mixture, input_data, prefix)
        type(Long_Interactions_Wrapper), intent(out) :: long_interactions
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        logical :: are_dipolar(size(mixture%components))

        call long_interactions_set(are_dipolar, mixture%components)
        call long_interactions_check(any(are_dipolar), environment%reciprocal_lattice, &
            environment%permittivity)
        call long_interactions_create(long_interactions%alpha, environment%periodic_box, &
            any(are_dipolar), input_data, prefix)

        call ewald_real_create(long_interactions%real_visitor, environment%periodic_box, &
            any(are_dipolar))
        call ewald_real_create(long_interactions%real_pairs, environment%permittivity, mixture%&
            components_min_distances, are_dipolar, long_interactions%alpha, input_data, &
            prefix//"Real.")
        call ewald_real_create(long_interactions%real_components, environment%periodic_box, &
            mixture%components, long_interactions%real_pairs)

        call ewald_reci_create(long_interactions%reci_weight, environment, long_interactions%&
            alpha, any(are_dipolar))
        call ewald_reci_create(long_interactions%reci_structures, environment, mixture%components, &
            are_dipolar, long_interactions%reci_weight)
    end subroutine create_all

    subroutine destroy_all(long_interactions)
        type(Long_Interactions_Wrapper), intent(inout) :: long_interactions

        call ewald_reci_destroy(long_interactions%reci_structures)
        call ewald_reci_destroy(long_interactions%reci_weight)

        call ewald_real_destroy(long_interactions%real_components)
        call ewald_real_destroy(long_interactions%real_pairs)
        call ewald_real_destroy(long_interactions%real_visitor)

        call long_interactions_destroy(long_interactions%alpha)
    end subroutine destroy_all

    subroutine create_alpha(alpha, periodic_box, dipoles_exist, input_data, prefix)
        class(Abstract_Ewald_Convergence_Parameter), allocatable, intent(out) :: alpha
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
            allocate(Concrete_Ewald_Convergence_Parameter :: alpha)
        else
            allocate(Null_Ewald_Convergence_Parameter :: alpha)
        end if
        call alpha%construct(periodic_box, alpha_x_box)
    end subroutine create_alpha

    subroutine destroy_alpha(alpha)
        class(Abstract_Ewald_Convergence_Parameter), allocatable, intent(inout) :: alpha

        if (allocated(alpha)) then
            call alpha%destroy()
            deallocate(alpha)
        end if
    end subroutine destroy_alpha

    subroutine set_are_dipolar(are_dipolar, components)
        logical, intent(out) :: are_dipolar(:)
        type(Component_Wrapper), intent(in) :: components(:)

        integer :: i_component

        do i_component = 1, size(are_dipolar)
            are_dipolar(i_component) = component_is_dipolar(components(i_component)%dipolar_moments)
        end do
    end subroutine set_are_dipolar

    subroutine check_consistency(dipoles_exist, reciprocal_lattice, permittivity)
        logical, intent(in) :: dipoles_exist
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Permittivity), intent(in) :: permittivity

        if (dipoles_exist) then
            if (.not.use_reciprocal_lattice(reciprocal_lattice)) then
                call warning_continue("check_consistency: dipoles exist but reciprocal_lattice "//&
                    "unused.")
            end if
            if (.not.use_permittivity(permittivity)) then
                call warning_continue("check_consistency: dipoles exist but permittivity unused.")
            end if
        end if
    end subroutine check_consistency

end module procedures_long_interactions_factory
