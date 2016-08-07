module classes_box_volume_change

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_visit_condition, only: visit_condition_in_range => in_range, &
    visit_condition_lower => lower, visit_condition_unconditional => unconditional
use procedures_dipoles_field_interaction, only: dipoles_field_visit_component => visit_component
use classes_changed_box_size, only: Abstract_Changed_Box_Size
use types_reals_line, only: Reals_Line
use types_temporary_observables, only: Concrete_Single_Energies, Concrete_Energies
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use procedures_triangle_observables, only: triangle_observables_sum
use procedures_generating_observables_factory, only:generating_observables_create => create
use classes_metropolis_algorithm, only: Abstract_Metropolis_Algorithm

implicit none

private

    type, extends(Abstract_Metropolis_Algorithm), abstract, public :: Abstract_Box_Volume_Change
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        class(Abstract_Changed_Box_Size), pointer :: changed_box_size => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: get_num_choices => Abstract_get_num_choices
        procedure :: try => Abstract_try
        procedure :: test_metropolis => Abstract_test_metropolis
        procedure :: visit_walls => Abstrac_visit_walls
        procedure :: visit_fields => Abstract_visit_fields
        procedure :: visit_short => Abstract_visit_short
    end type Abstract_Box_Volume_Change

contains

!implementation Abstract_Metropolis_Algorithm

    subroutine Abstract_construct(this, environment, mixture, short_interactions, changed_box_size)
        class(Abstract_Box_Volume_Change), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        class(Abstract_Changed_Box_Size), target, intent(in) :: changed_box_size

        this%environment => environment
        this%mixture => mixture
        this%short_interactions => short_interactions
        this%changed_box_size => changed_box_size
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Box_Volume_Change), intent(inout) :: this

        this%changed_box_size => null()
        this%short_interactions => null()
        this%mixture => null()
        this%environment => null()
    end subroutine Abstract_destroy

    pure integer function Abstract_get_num_choices(this) result(num_choices)
        class(Abstract_Box_Volume_Change), intent(in) :: this

        num_choices = 1
    end function Abstract_get_num_choices

    subroutine Abstract_try(this, observables)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Energies) :: deltas

        observables%volume_change_counter%num_hits = observables%volume_change_counter%num_hits + 1
        call generating_observables_create(deltas%walls_energies, size(observables%walls_energies))
        call generating_observables_create(deltas%field_energies, size(observables%field_energies))
        call this%test_metropolis(success, deltas)
    end subroutine Abstract_try

    subroutine Abstract_test_metropolis(this, success, deltas)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Energies), intent(inout) :: deltas

        real(DP) :: delta_energy
        type(Concrete_Energies) :: new_energies
        real(DP) :: box_size_ratio(num_dimensions)
        integer :: i_component
        logical :: overlap

        success = .false.
        box_size_ratio = this%changed_box_size%get_ratio()
        do i_component = 1, size(this%mixture%components)
            call this%mixture%components(i_component)%positions%rescale_all(box_size_ratio)
        end do

        call generating_observables_create(new_energies%walls_energies, size(deltas%walls_energies))
        call this%visit_walls(overlap, new_energies%walls_energies)
        if (overlap) return
        call generating_observables_create(new_energies%field_energies, size(deltas%field_energies))
        call this%visit_fields(new_energies%field_energies)
        call generating_observables_create(new_energies%short_energies, size(deltas%short_energies))
        call this%visit_short(overlap, new_energies%short_energies)
        if (overlap) return
        !dipolar interactions

        delta_energy = sum(deltas%walls_energies + deltas%field_energies) + &
            triangle_observables_sum(deltas%short_energies)
    end subroutine Abstract_test_metropolis

    subroutine Abstrac_visit_walls(this, overlap, new_energies)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: new_energies(:)

        integer :: i_component

        do i_component = 1, size(new_energies)
            call this%short_interactions%walls_visitor%visit(overlap, new_energies(i_component), &
                this%mixture%components(i_component)%positions, this%short_interactions%&
                wall_pairs(i_component)%potential)
            if (overlap) return
        end do
    end subroutine Abstrac_visit_walls

    subroutine Abstract_visit_fields(this, new_energies)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        real(DP), intent(out) :: new_energies(:)

        integer :: i_component

        do i_component = 1, size(new_energies)
            new_energies(i_component) = dipoles_field_visit_component(this%environment%&
                external_field, this%mixture%components(i_component)%positions, this%mixture%&
                components(i_component)%dipolar_moments)
        end do
    end subroutine Abstract_visit_fields

    subroutine Abstract_visit_short(this, overlap, new_energies)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        logical, intent(out) :: overlap
        type(Reals_Line), intent(inout) :: new_energies(:)

        real(DP) :: energy_ij, energy_j
        integer :: j_component, i_component, i_particle, i_exclude
        logical :: same_component
        type(Concrete_Temporary_Particle) :: particle
        procedure(visit_condition_in_range), pointer :: in_range => null()

        do j_component = 1, size(this%short_interactions%visitable_cells, 2)
            do i_component = 1, j_component
                same_component = i_component == j_component
                if (same_component) then
                    in_range => visit_condition_lower
                else
                    in_range => visit_condition_unconditional
                end if
                energy_ij = 0._DP
                do i_particle = 1, this%mixture%components(j_component)%positions%get_num()
                    particle%i = i_particle
                    particle%position = this%mixture%components(j_component)%positions%&
                        get(particle%i)
                    i_exclude = merge(particle%i, 0, same_component)
                    call this%short_interactions%visitable_cells(i_component, j_component)%&
                        visit_energy(overlap, energy_j, particle, in_range, i_exclude)
                    if (overlap) return
                    energy_ij = energy_ij + energy_j
                end do
                new_energies(j_component)%line(i_component) = energy_ij
            end do
        end do
    end subroutine Abstract_visit_short

!end implementation Abstract_Metropolis_Algorithm

end module classes_box_volume_change
