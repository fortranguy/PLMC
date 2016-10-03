module procedures_dipolar_interactions_factory

use json_module, only: json_file
use procedures_errors, only: warning_continue
use procedures_boxes_factory, only: boxes_create => create, boxes_destroy => destroy
use classes_permittivity, only: Abstract_Permittivity
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_inquirers, only: use_permittivity, use_reciprocal_lattice
use types_mixture_wrapper, only: Mixture_Wrapper
use procedures_mixture_total_moments_factory, only: set_are_dipolar
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use procedures_des_convergence_parameter_factory, only: des_convergence_parameter_create => create,&
    des_convergence_parameter_destroy => destroy
use procedures_des_real_factory, only: des_real_create => create, des_real_destroy => destroy
use procedures_des_reci_factory, only: des_reci_create => create, des_reci_destroy => destroy
use procedures_des_self_factory, only: des_self_create => create, des_self_destroy => destroy
use procedures_des_surf_mixture_factory, only: des_surf_mixture_create => create, &
    des_surf_mixture_destroy => destroy
use procedures_dlc_factory, only: dlc_create => create, dlc_destroy => destroy

implicit none

private
public :: dipolar_interactions_create, dipolar_interactions_destroy

contains

    subroutine dipolar_interactions_create(dipolar_interactions_dynamic, &
        dipolar_interactions_static, environment, mixture, generating_data, prefix)
        type(Dipolar_Interactions_Dynamic_Wrapper), intent(out) :: dipolar_interactions_dynamic
        type(Dipolar_Interactions_Static_Wrapper), allocatable, intent(out) :: &
            dipolar_interactions_static(:)
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: i_box
        logical :: are_dipolar(size(mixture%gemc_components, 1), size(mixture%gemc_components, 2))

        do i_box = 1, size(are_dipolar, 2)
            call set_are_dipolar(are_dipolar(:, i_box), mixture%gemc_components(:, i_box))
            call check_consistency(environment%reciprocal_lattices(i_box), environment%&
                permittivity, any(are_dipolar(:, i_box)))
        end do

        call des_convergence_parameter_create(dipolar_interactions_dynamic%alpha, any(are_dipolar),&
            generating_data, prefix)

        allocate(dipolar_interactions_static(size(mixture%gemc_components, 2)))

        do i_box = 1, size(dipolar_interactions_static)
            call boxes_create(dipolar_interactions_static(i_box)%box_size_memento_real, &
                environment%periodic_boxes(i_box), environment%beta_pressure, &
                any(are_dipolar(:, i_box)))
            call des_real_create(dipolar_interactions_static(i_box)%real_pair, &
                dipolar_interactions_static(i_box)%box_size_memento_real, environment%permittivity,&
                mixture%components_min_distances, any(are_dipolar(:, i_box)), &
                dipolar_interactions_dynamic%alpha, generating_data, prefix//"Real.")
        end do
        call des_real_create(dipolar_interactions_dynamic%gemc_real_components, environment%&
            periodic_boxes, mixture%gemc_components, are_dipolar, dipolar_interactions_static)

        do i_box = 1, size(dipolar_interactions_static)
            allocate(dipolar_interactions_static(i_box)%box_size_memento_reci, &
                source=dipolar_interactions_static(i_box)%box_size_memento_real)
            call dipolar_interactions_static(i_box)%box_size_memento_reci%target(environment%&
                periodic_boxes(i_box))
            call des_reci_create(dipolar_interactions_static(i_box)%reci_weight, &
                dipolar_interactions_static(i_box)%box_size_memento_reci, environment%&
                reciprocal_lattices(i_box), environment%permittivity, any(are_dipolar(:, i_box)), &
                dipolar_interactions_dynamic%alpha)
            call des_reci_create(dipolar_interactions_static(i_box)%reci_structure, environment%&
                periodic_boxes(i_box), dipolar_interactions_static(i_box)%box_size_memento_reci, &
                environment%reciprocal_lattices(i_box), mixture%gemc_components(:, i_box), &
                are_dipolar(:, i_box))
        end do
        call des_reci_create(dipolar_interactions_dynamic%reci_visitors, environment%&
                periodic_boxes, environment%reciprocal_lattices, dipolar_interactions_static)

        call des_self_create(dipolar_interactions_dynamic%gemc_self_components, environment%&
            periodic_boxes, environment%permittivity, mixture%gemc_components, are_dipolar, &
            dipolar_interactions_dynamic%alpha)

        call des_surf_mixture_create(dipolar_interactions_dynamic%gemc_surf_mixture, environment%&
            periodic_boxes, environment%permittivity, mixture%total_moments)

        do i_box = 1, size(dipolar_interactions_static)
            call dlc_create(dipolar_interactions_static(i_box)%dlc_weight, environment%&
                periodic_boxes(i_box), environment%reciprocal_lattices(i_box), environment%&
                permittivity, any(are_dipolar(:, i_box)))
            call dlc_create(dipolar_interactions_static(i_box)%dlc_structures, environment%&
                periodic_boxes(i_box), environment%reciprocal_lattices(i_box), mixture%&
                gemc_components(:, i_box), are_dipolar(:, i_box))
        end do
        call dlc_create(dipolar_interactions_dynamic%dlc_visitors, environment%&
                periodic_boxes, environment%reciprocal_lattices, dipolar_interactions_static)
    end subroutine dipolar_interactions_create

    !> @todo
    !> if (allocated(dipolar_interactions_static)): before?
    subroutine dipolar_interactions_destroy(dipolar_interactions_dynamic, &
        dipolar_interactions_static)
        type(Dipolar_Interactions_Dynamic_Wrapper), intent(inout) :: dipolar_interactions_dynamic
        type(Dipolar_Interactions_Static_Wrapper), allocatable, intent(inout) :: &
            dipolar_interactions_static(:)

        integer :: i_box

        call dlc_destroy(dipolar_interactions_dynamic%dlc_visitors)
        do i_box = size(dipolar_interactions_static), 1, -1
            call dlc_destroy(dipolar_interactions_static(i_box)%dlc_structures)
            call dlc_destroy(dipolar_interactions_static(i_box)%dlc_weight)
        end do

        call des_surf_mixture_destroy(dipolar_interactions_dynamic%gemc_surf_mixture)

        call des_self_destroy(dipolar_interactions_dynamic%gemc_self_components)

        call des_reci_destroy(dipolar_interactions_dynamic%reci_visitors)
        do i_box = size(dipolar_interactions_static), 1, -1
            call des_reci_destroy(dipolar_interactions_static(i_box)%reci_structure)
            call des_reci_destroy(dipolar_interactions_static(i_box)%reci_weight)
            call boxes_destroy(dipolar_interactions_static(i_box)%box_size_memento_reci)
        end do

        call des_real_destroy(dipolar_interactions_dynamic%gemc_real_components)
        do i_box = size(dipolar_interactions_static), 1, -1
            call des_real_destroy(dipolar_interactions_static(i_box)%real_pair)
            call boxes_destroy(dipolar_interactions_static(i_box)%box_size_memento_real)
        end do

        if (allocated(dipolar_interactions_static)) deallocate(dipolar_interactions_static)

        call des_convergence_parameter_destroy(dipolar_interactions_dynamic%alpha)
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
