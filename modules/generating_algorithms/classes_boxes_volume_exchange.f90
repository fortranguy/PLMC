module classes_boxes_volume_exchange

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use types_logical_line, only: Logical_Triangle
use procedures_logical_factory, only: logical_create => create
use classes_hetero_couples, only: Abstract_Hetero_Couples
use procedures_hetero_couples_factory, only: hetero_couples_destroy => destroy
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only: tower_sampler_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_factory, only: mixture_rescale_positions => rescale_positions
use types_cells_wrapper, only: Cells_Wrapper
use procedures_cells_memento, only: cells_memento_save => save, cells_memento_restore => restore
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_short_interactions_resetter, only: short_interactions_reset => reset
use procedures_short_interactions_visitor, only: short_interactions_visit => visit
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use procedures_dipolar_interactions_factory, only: dipolar_interactions_destroy => destroy
use procedures_dipolar_interactions_visitor, only: dipolar_interactions_visit => visit
use classes_dipolar_interactions_facade, only: Abstract_Dipolar_Interactions_Facade
use classes_exchanged_boxes_size, only: Exchanged_Boxes_Size_Line
use procedures_triangle_observables, only: triangle_observables_diff, triangle_observables_sum
use types_observables_energies, only: Concrete_Observables_Energies
use procedures_observables_energies_factory, only: observables_energies_create => create, &
    observables_energies_set => set
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use classes_generating_algorithm, only: Abstract_Generating_Algorithm
use procedures_selectors_resetters, only: selectors_reset => reset
use procedures_metropolis_algorithm, only: metropolis_algorithm

implicit none

private

    type, extends(Abstract_Generating_Algorithm), public :: Boxes_Volume_Exchange
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Component_Wrapper), pointer :: components(:, :) => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        class(Abstract_Dipolar_Interactions_Facade), pointer :: dipolar_interactions_facades(:) => &
            null()
        type(Exchanged_Boxes_Size_Line), pointer :: exchanged_boxes_size(:) => null()
        logical, allocatable :: have_positions(:, :)
        class(Abstract_Hetero_Couples), allocatable :: couples
        class(Abstract_Tower_Sampler), allocatable :: selector
    contains
        procedure :: construct => Concrete_construct
        procedure :: destroy => Concrete_destroy
        procedure :: reset_selectors => Concrete_reset_selectors
        procedure :: get_num_choices => Concrete_get_num_choices
        procedure :: try => Concrete_try
        procedure, private :: metropolis_algorithm => Concrete_metropolis_algorithm
        procedure, private :: acceptation_probability => Concrete_acceptation_probability
    end type Boxes_Volume_Exchange

contains

    subroutine Concrete_construct(this, environment, components, short_interactions, &
        dipolar_interactions_facades, exchanged_boxes_size, have_positions, couples, selector)
        class(Boxes_Volume_Exchange), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:, :)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        class(Abstract_Dipolar_Interactions_Facade), target, intent(in) :: &
            dipolar_interactions_facades(:)
        type(Exchanged_Boxes_Size_Line), target, intent(in) :: exchanged_boxes_size(:)
        logical, intent(in) :: have_positions(:, :)
        class(Abstract_Hetero_Couples), intent(in) :: couples
        class(Abstract_Tower_Sampler), intent(in) :: selector

        this%environment => environment
        this%components => components
        this%short_interactions => short_interactions
        this%dipolar_interactions_facades => dipolar_interactions_facades
        this%exchanged_boxes_size => exchanged_boxes_size
        allocate(this%have_positions, source=have_positions)
        allocate(this%couples, source=couples)
        allocate(this%selector, source=selector)
    end subroutine Concrete_construct

    subroutine Concrete_destroy(this)
        class(Boxes_Volume_Exchange), intent(inout) :: this

        call tower_sampler_destroy(this%selector)
        call hetero_couples_destroy(this%couples)
        if (allocated(this%have_positions)) deallocate(this%have_positions)
        this%exchanged_boxes_size => null()
        this%dipolar_interactions_facades => null()
        this%short_interactions => null()
        this%components => null()
        this%environment => null()
    end subroutine Concrete_destroy

    subroutine Concrete_reset_selectors(this)
        class(Boxes_Volume_Exchange), intent(inout) :: this

        call selectors_reset(this%selector, this%couples, this%exchanged_boxes_size, this%&
            components, this%have_positions)
    end subroutine Concrete_reset_selectors

    pure integer function Concrete_get_num_choices(this) result(num_choices)
        class(Boxes_Volume_Exchange), intent(in) :: this

        num_choices = this%selector%get_num_choices()
    end function Concrete_get_num_choices

    !> @note About [[procedures_dipolar_interactions_factory:destroy_static]], cf.
    !> [[classes_box_volume_change:Concrete_try]]
    subroutine Concrete_try(this, observables)
        class(Boxes_Volume_Exchange), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Observables_Energies) :: new_energies(2)
        real(DP), dimension(num_dimensions, size(new_energies)) :: new_boxes_size, boxes_size, &
            boxes_size_ratio
        integer :: ij_boxes(size(new_energies)), i_box, j_box, i_partner
        type(Cells_Wrapper) :: cells(size(new_energies))
        type(Logical_Triangle)  :: only_resized_triangle(size(new_energies))
        type(Dipolar_Interactions_Static_Wrapper) :: dipolar_interactions_static(size(new_energies))
        logical :: reset_real_pair(size(new_energies))

        ij_boxes = this%couples%get(this%selector%get())
        j_box = maxval(ij_boxes)
        i_box = minval(ij_boxes)

        observables%volumes_exchange_counter(j_box)%line(i_box)%num_hits = &
            observables%volumes_exchange_counter(j_box)%line(i_box)%num_hits + 1

        do i_partner = 1, size(new_boxes_size, 2)
            boxes_size(:, i_partner) = this%environment%periodic_boxes(ij_boxes(i_partner))%&
                get_size()
        end do
        boxes_size_ratio = this%exchanged_boxes_size(j_box)%line(i_box)%exchanged%&
            get_ratios(boxes_size(:, 1) / boxes_size(:, 2))
        new_boxes_size = boxes_size * boxes_size_ratio

        do i_partner = 1, size(dipolar_interactions_static)
            call this%dipolar_interactions_facades(ij_boxes(i_partner))%&
                save(dipolar_interactions_static(i_partner), reset_real_pair(i_partner), &
                    new_boxes_size(:, i_partner))
        end do

        do i_partner = 1, size(new_boxes_size, 2)
            call this%environment%periodic_boxes(ij_boxes(i_partner))%&
                set(new_boxes_size(:, i_partner))
            call mixture_rescale_positions(this%components(:, ij_boxes(i_partner)), &
                boxes_size_ratio(:, i_partner))
            call logical_create(only_resized_triangle(i_partner)%triangle, size(this%components, 1))
            call cells_memento_save(cells(i_partner), only_resized_triangle(i_partner)%triangle, &
                this%short_interactions%visitable_cells_memento, this%short_interactions%&
                cells(ij_boxes(i_partner)))
            call short_interactions_reset(this%short_interactions%cells(ij_boxes(i_partner))%&
                neighbour_cells, only_resized_triangle(i_partner)%triangle, this%&
                short_interactions%cells(ij_boxes(i_partner))%visitable_cells)
            call this%dipolar_interactions_facades(ij_boxes(i_partner))%&
                reset(reset_real_pair(i_partner))
            call observables_energies_create(new_energies(i_partner), size(this%components, 1))
        end do

        call this%metropolis_algorithm(success, new_energies, ij_boxes, &
            product(boxes_size_ratio, 1), observables%energies)

        if (success) then
            do i_partner = 1, size(ij_boxes)
                observables%accessible_domains_size(:, ij_boxes(i_partner)) = this%environment%&
                    accessible_domains(ij_boxes(i_partner))%get_size()
                call observables_energies_set(observables%energies(ij_boxes(i_partner)), &
                    new_energies(i_partner))
            end do
            observables%volumes_exchange_counter(j_box)%line(i_box)%num_successes = &
                observables%volumes_exchange_counter(j_box)%line(i_box)%num_successes + 1
        else
            do i_partner = 1, size(ij_boxes)
                call this%environment%periodic_boxes(ij_boxes(i_partner))%&
                    set(boxes_size(:, i_partner))
                call mixture_rescale_positions(this%components(:, ij_boxes(i_partner)), &
                    1._DP / boxes_size_ratio(:, i_partner))
                call cells_memento_restore(this%short_interactions%cells(ij_boxes(i_partner)), &
                    only_resized_triangle(i_partner)%triangle, this%short_interactions%&
                    visitable_cells_memento, cells(i_partner))
                call this%dipolar_interactions_facades(ij_boxes(i_partner))%&
                    restore(dipolar_interactions_static(i_partner), reset_real_pair(i_partner))
            end do
        end if

        do i_partner = size(dipolar_interactions_static), 1, -1
            call dipolar_interactions_destroy(dipolar_interactions_static(i_partner))
        end do
    end subroutine Concrete_try

    subroutine Concrete_metropolis_algorithm(this, success, new_energies, ij_boxes, &
        boxes_volume_ratio, energies)
        class(Boxes_Volume_Exchange), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Observables_Energies), intent(inout) :: new_energies(:)
        integer, intent(in) :: ij_boxes(:)
        real(DP), intent(in) :: boxes_volume_ratio(:)
        type(Concrete_Observables_Energies), intent(in) :: energies(:)

        real(DP) :: delta_energies(size(new_energies))
        type(Concrete_Observables_Energies) :: deltas(size(new_energies))
        logical :: overlap
        integer :: i_partner

        success = .false.
        do i_partner = 1, size(deltas)
            call observables_energies_create(deltas(i_partner), size(this%components, 1))
        end do

        do i_partner = 1, size(deltas)
            call short_interactions_visit(overlap, new_energies(i_partner)%walls_energies, this%&
                components(:, ij_boxes(i_partner)), this%short_interactions%&
                walls_visitors(ij_boxes(i_partner)), this%short_interactions%wall_pairs)
            if (overlap) return
            deltas(i_partner)%walls_energies = new_energies(i_partner)%walls_energies - &
                energies(ij_boxes(i_partner))%walls_energies
        end do

        do i_partner = 1, size(deltas)
            call short_interactions_visit(overlap, new_energies(i_partner)%short_energies, this%&
                components(:, ij_boxes(i_partner)), this%short_interactions%&
                cells(ij_boxes(i_partner))%visitable_cells)
            if (overlap) return
            call triangle_observables_diff(deltas(i_partner)%short_energies, &
                new_energies(i_partner)%short_energies, energies(ij_boxes(i_partner))%&
                short_energies)
        end do

        do i_partner = 1, size(deltas)
            call dipolar_interactions_visit(new_energies(i_partner)%field_energies, this%&
                environment%external_fields(ij_boxes(i_partner)), this%&
                components(:, ij_boxes(i_partner)))
            deltas(i_partner)%field_energies = new_energies(i_partner)%field_energies - &
                energies(ij_boxes(i_partner))%field_energies
            call this%dipolar_interactions_facades(ij_boxes(i_partner))%&
                visit(new_energies(i_partner)%dipolar_energies, new_energies(i_partner)%&
                    dipolar_shared_energy, boxes_volume_ratio(i_partner), &
                    energies(ij_boxes(i_partner))%dipolar_energies, energies(ij_boxes(i_partner))%&
                    dipolar_shared_energy)
            call triangle_observables_diff(deltas(i_partner)%dipolar_energies, &
                new_energies(i_partner)%dipolar_energies, energies(ij_boxes(i_partner))%&
                dipolar_energies)
            deltas(i_partner)%dipolar_shared_energy = new_energies(i_partner)%&
                dipolar_shared_energy - energies(ij_boxes(i_partner))%dipolar_shared_energy
        end do

        do i_partner = 1, size(delta_energies)
            delta_energies(i_partner) = &
                sum(deltas(i_partner)%walls_energies + deltas(i_partner)%field_energies) + &
                triangle_observables_sum(deltas(i_partner)%short_energies) + &
                triangle_observables_sum(deltas(i_partner)%dipolar_energies) + &
                deltas(i_partner)%dipolar_shared_energy
        end do
        success = metropolis_algorithm(this%acceptation_probability(ij_boxes, boxes_volume_ratio, &
            sum(delta_energies)))
    end subroutine Concrete_metropolis_algorithm

    !> \[
    !>      P[(V_{\boldsymbol{I}}, V_{\boldsymbol{J}}) \to
    !>          (V_{\boldsymbol{I}}^\prime, V_{\boldsymbol{J}}^\prime)] = \min \left[ 1,
    !>              \left( \frac{V_{\boldsymbol{I}}^\prime}
    !>                  {V_{\boldsymbol{I}}} \right)^{N_{\boldsymbol{I}} + 1}
    !>              \left( \frac{V_{\boldsymbol{J}}^\prime}
    !>                  {V_{\boldsymbol{J}}} \right)^{N_{\boldsymbol{J}} + 1}
    !>              e^{-\beta \Delta U_{(V_{\boldsymbol{I}}, V_{\boldsymbol{J}}) \to
    !>                  (V_{\boldsymbol{I}}^\prime, V_{\boldsymbol{J}}^\prime)}}
    !>          \right]
    !> \]
    pure real(DP) function Concrete_acceptation_probability(this, ij_boxes, boxes_volume_ratio, &
        delta_energy) result(probability)
        class(Boxes_Volume_Exchange), intent(in) :: this
        integer, intent(in) :: ij_boxes(:)
        real(DP), intent(in) :: boxes_volume_ratio(:)
        real(DP), intent(in) :: delta_energy

        integer :: i_component, nums_particles(size(ij_boxes)), i_partner

        do i_partner = 1, size(nums_particles)
            nums_particles(i_partner) = 0
            do i_component = 1, size(this%components, 1)
                nums_particles(i_partner) = nums_particles(i_partner) + this%&
                    components(i_component, ij_boxes(i_partner))%num_particles%get()
            end do
        end do

        probability = &
            boxes_volume_ratio(1)**(nums_particles(1) + 1) * &
            boxes_volume_ratio(2)**(nums_particles(2) + 1) * &
            exp(-delta_energy / this%environment%temperature%get())
        probability = min(1._DP, probability)
    end function Concrete_acceptation_probability

end module classes_boxes_volume_exchange
