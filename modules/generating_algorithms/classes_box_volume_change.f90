module classes_box_volume_change

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use types_logical_line, only: Concrete_Logical_Line
use procedures_logical_factory, only: logical_create => create
use procedures_random_number, only: random_integer
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only: tower_sampler_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use classes_neighbour_cells, only: Neighbour_Cells_Line
use classes_visitable_cells, only: Abstract_Visitable_Cells
use procedures_cells_factory, only: cells_destroy => destroy, cells_allocate_triangle => &
    allocate_triangle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_short_interactions_resetter, only: short_interactions_reset => reset
use procedures_short_interactions_visitor, only: short_interactions_visit => visit, &
    short_interactions_visit_cells => visit_cells
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use procedures_dipolar_interactions_visitor, only: dipolar_interactions_visit => visit
use classes_dipolar_interactions_facade, only: Abstract_Dipolar_Interactions_Facade
use classes_changed_box_size, only: Changed_Box_Size_Line
use procedures_triangle_observables, only: operator(-)
use types_observables_energies, only: Concrete_Observables_Energies
use procedures_observables_energies_factory, only: observables_energies_set => set
use procedures_triangle_observables, only: triangle_observables_sum
use procedures_observables_energies_factory, only: observables_energies_create => create
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use classes_generating_algorithm, only: Abstract_Generating_Algorithm
use procedures_selectors_resetters, only: selectors_reset => reset
use procedures_metropolis_algorithm, only: metropolis_algorithm

implicit none

private

    type, extends(Abstract_Generating_Algorithm), public :: Box_Volume_Change
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        class(Abstract_Dipolar_Interactions_Facade), pointer :: dipolar_interactions_facades(:) => &
            null()
        type(Changed_Box_Size_Line), pointer :: changed_boxes_size(:) => null()
        logical, allocatable :: have_positions(:, :)
        class(Abstract_Tower_Sampler), allocatable :: selectors(:)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset_selectors => Abstract_reset_selectors
        procedure :: get_num_choices => Abstract_get_num_choices
        procedure :: try => Abstract_try
        procedure, private :: metropolis_algorithm => Abstract_metropolis_algorithm
        procedure, private :: acceptation_probability => Abstract_acceptation_probability
        procedure, private :: rescale_positions => Abstract_rescale_positions
        procedure, private :: save_cells => Abstract_save_cells
        procedure, private :: restore_cells => Abstract_restore_cells
    end type Box_Volume_Change

contains

    subroutine Abstract_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_facades, changed_boxes_size, have_positions, selectors)
        class(Box_Volume_Change), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        class(Abstract_Dipolar_Interactions_Facade), target, intent(in) :: &
            dipolar_interactions_facades(:)
        type(Changed_Box_Size_Line), target, intent(in) :: changed_boxes_size(:)
        logical, intent(in) :: have_positions(:, :)
        class(Abstract_Tower_Sampler), intent(in) :: selectors(:)

        this%environment => environment
        this%mixture => mixture
        this%short_interactions => short_interactions
        this%dipolar_interactions_facades => dipolar_interactions_facades
        this%changed_boxes_size => changed_boxes_size
        allocate(this%have_positions, source=have_positions)
        allocate(this%selectors, source=selectors)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Box_Volume_Change), intent(inout) :: this

        call tower_sampler_destroy(this%selectors)
        if (allocated(this%have_positions)) deallocate(this%have_positions)
        this%changed_boxes_size => null()
        this%dipolar_interactions_facades => null()
        this%short_interactions => null()
        this%mixture => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_reset_selectors(this)
        class(Box_Volume_Change), intent(inout) :: this

        call selectors_reset(this%selectors, this%changed_boxes_size, this%mixture%components,&
            this%have_positions)
    end subroutine Abstract_reset_selectors

    pure integer function Abstract_get_num_choices(this) result(num_choices)
        class(Box_Volume_Change), intent(in) :: this

        integer :: i_box

        num_choices = 0
        do i_box = 1, size(this%selectors)
            num_choices = num_choices + this%selectors(i_box)%get_num_choices()
        end do
    end function Abstract_get_num_choices

    subroutine Abstract_try(this, observables)
        class(Box_Volume_Change), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Observables_Energies) :: new_energies
        real(DP), dimension(num_dimensions) :: new_box_size, box_size, box_size_ratio
        integer :: i_box
        type(Neighbour_Cells_Line), allocatable :: neighbour_cells(:)
        type(Concrete_Logical_Line), allocatable :: only_resized_triangle(:)
        class(Abstract_Visitable_Cells), allocatable :: visitable_cells(:, :)
        type(Dipolar_Interactions_Static_Wrapper) :: dipolar_interactions_static
        logical :: reset_real_pair

        i_box = random_integer(size(this%environment%periodic_boxes))
        observables%volumes_change_counter(i_box)%line(i_box)%num_hits = observables%&
            volumes_change_counter(i_box)%line(i_box)%num_hits + 1
        box_size = this%environment%periodic_boxes(i_box)%get_size()
        box_size_ratio = this%changed_boxes_size(i_box)%line(i_box)%changed%get_ratio()
        new_box_size = box_size * box_size_ratio
        call this%dipolar_interactions_facades(i_box)%save(dipolar_interactions_static, reset_real_pair, new_box_size)

        call this%environment%periodic_boxes(i_box)%set(new_box_size)
        call this%rescale_positions(i_box, box_size_ratio)
        call logical_create(only_resized_triangle, size(this%mixture%components, 1))
        call this%save_cells(neighbour_cells, only_resized_triangle, visitable_cells, i_box)
        call short_interactions_reset(this%short_interactions%cells(i_box)%neighbour_cells, &
            only_resized_triangle, this%short_interactions%cells(i_box)%visitable_cells)
        call this%dipolar_interactions_facades(i_box)%reset()
        call observables_energies_create(new_energies, size(this%mixture%components, 1))
        call this%metropolis_algorithm(success, new_energies, i_box, box_size_ratio, observables%energies(i_box))

        if (success) then
            observables%accessible_domains_size(:, i_box) = this%environment%&
                accessible_domains(i_box)%get_size()
            call observables_energies_set(observables%energies(i_box), new_energies)
            observables%volumes_change_counter(i_box)%line(i_box)%num_successes = observables%&
                volumes_change_counter(i_box)%line(i_box)%num_successes + 1
        else
            call this%environment%periodic_boxes(i_box)%set(box_size)
            call this%rescale_positions(i_box, 1._DP / box_size_ratio)
            call this%restore_cells(neighbour_cells, only_resized_triangle, visitable_cells, i_box)
            call this%dipolar_interactions_facades(i_box)%restore(dipolar_interactions_static, reset_real_pair)
        end if
    end subroutine Abstract_try

    subroutine Abstract_metropolis_algorithm(this, success, new_energies, i_box, box_size_ratio, energies)
        class(Box_Volume_Change), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Observables_Energies), intent(inout) :: new_energies
        integer, intent(in) :: i_box
        real(DP), intent(in) :: box_size_ratio(:)
        type(Concrete_Observables_Energies), intent(in) :: energies

        real(DP) :: delta_energy
        type(Concrete_Observables_Energies) :: deltas

        logical :: overlap

        success = .false.
        call observables_energies_create(deltas, size(this%mixture%components, 1))
        call short_interactions_visit(overlap, new_energies%walls_energies, this%mixture%&
            components(:, i_box), this%short_interactions%walls_visitors(i_box), this%&
            short_interactions%wall_pairs)
        if (overlap) return
        deltas%walls_energies = new_energies%walls_energies - energies%walls_energies
        call short_interactions_visit_cells(overlap, new_energies%short_energies, this%mixture%&
            components(:, i_box), this%short_interactions%cells(i_box)%visitable_cells)
        deltas%short_energies = new_energies%short_energies - energies%short_energies
        if (overlap) return
        call dipolar_interactions_visit(new_energies%field_energies, this%environment%&
            external_fields(i_box), this%mixture%components(:, i_box))
        deltas%field_energies = new_energies%field_energies - energies%field_energies
        call this%dipolar_interactions_facades(i_box)%visit(new_energies%dipolar_energies, new_energies%&
            dipolar_shared_energy, product(box_size_ratio), energies%dipolar_energies, energies%&
            dipolar_shared_energy)
        deltas%dipolar_energies = new_energies%dipolar_energies - energies%dipolar_energies
        deltas%dipolar_shared_energy = new_energies%dipolar_shared_energy - energies%&
            dipolar_shared_energy

        delta_energy = sum(deltas%walls_energies + deltas%field_energies) + &
            triangle_observables_sum(deltas%short_energies) + &
            triangle_observables_sum(deltas%dipolar_energies) + deltas%dipolar_shared_energy
        success = metropolis_algorithm(this%acceptation_probability(i_box, box_size_ratio, delta_energy))
    end subroutine Abstract_metropolis_algorithm

    !> \[
    !>      P[V \to V^\prime] = min \left( 1, e^{-\beta p V \left( \frac{V^\prime}{V} - 1 \right)}
    !>          \left( \frac{V^\prime}{V} \right)^{N+1}
    !>          e^{-\beta [U(\vec{s}^N, V^\prime) - U(\vec{s}^N, V)]} \right)
    !> \]
    pure real(DP) function Abstract_acceptation_probability(this, i_box, box_size_ratio, delta_energy) &
        result(probability)
        class(Box_Volume_Change), intent(in) :: this
        integer, intent(in) :: i_box
        real(DP), intent(in) :: box_size_ratio(:)
        real(DP), intent(in) :: delta_energy

        real(DP) :: volume_ratio
        integer :: i_component, num_particles

        num_particles = 0._DP
        do i_component = 1, size(this%mixture%components, 1)
            num_particles = num_particles + this%mixture%components(i_component, i_box)%num_particles%get()
        end do
        volume_ratio = product(box_size_ratio)

        probability = exp(-this%environment%beta_pressure%get() * &
            product(this%environment%accessible_domains(i_box)%get_size()) * (volume_ratio - 1._DP)) * &
            volume_ratio**(num_particles + 1) * &
            exp(-delta_energy / this%environment%temperature%get())
        probability = min(1._DP, probability)
    end function Abstract_acceptation_probability

    subroutine Abstract_rescale_positions(this, i_box, box_size_ratio)
        class(Box_Volume_Change), intent(in) :: this
        integer, intent(in) :: i_box
        real(DP), intent(in) :: box_size_ratio(:)

        integer :: i_component

        do i_component = 1, size(this%mixture%components, 1)
            call this%mixture%components(i_component, i_box)%positions%rescale_all(box_size_ratio)
        end do
    end subroutine Abstract_rescale_positions

    !> @note The methods [[classes_box_volume_change:Abstract_save_cells]] and
    !> [[classes_volume_change_method:Abstract_save_cells]]
    !> should be factorised. However, this would create a memory leak.
    !> @todo Send a bug report to gfortran.
    subroutine Abstract_save_cells(this, neighbour_cells, only_resized_triangle, visitable_cells, i_box)
        class(Box_Volume_Change), intent(in) :: this
        type(Neighbour_Cells_Line), allocatable, intent(out) :: neighbour_cells(:)
        type(Concrete_Logical_Line), intent(inout) ::only_resized_triangle(:)
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: visitable_cells(:, :)
        integer, intent(in) :: i_box

        integer :: i_component, j_component

        call cells_allocate_triangle(neighbour_cells, size(this%mixture%components, 1))
        do j_component = 1, size(neighbour_cells)
            do i_component = 1, size(neighbour_cells(j_component)%line)
                only_resized_triangle(j_component)%line(i_component) = this%short_interactions%cells(i_box)%&
                    neighbour_cells(j_component)%line(i_component)%cells%resize_only()
                if (.not. only_resized_triangle(j_component)%line(i_component)) then
                    allocate(neighbour_cells(j_component)%line(i_component)%cells, source=this%&
                        short_interactions%cells(i_box)%neighbour_cells(j_component)%line(i_component)%cells)
                end if
            end do
        end do
        call this%short_interactions%visitable_cells_memento%save(visitable_cells, this%&
            short_interactions%cells(i_box)%visitable_cells)
    end subroutine Abstract_save_cells

    !> @note The methods [[classes_box_volume_change:Abstract_restore_cells]] and
    !> [[classes_volume_change_method:Abstract_restore_cells]]
    !> should be factorised. cf. [[classes_box_volume_change:Abstract_save_cells]].
    subroutine Abstract_restore_cells(this, neighbour_cells, only_resized_triangle, visitable_cells, i_box)
        class(Box_Volume_Change), intent(in) :: this
        type(Neighbour_Cells_Line), allocatable, intent(inout) :: neighbour_cells(:)
        type(Concrete_Logical_Line), intent(in) ::only_resized_triangle(:)
        class(Abstract_Visitable_Cells), allocatable, intent(inout) :: visitable_cells(:, :)
        integer, intent(in) :: i_box

        integer :: i_component, j_component

        do j_component = 1, size(this%short_interactions%cells(i_box)%neighbour_cells)
            do i_component = 1, size(this%short_interactions%cells(i_box)%neighbour_cells(j_component)%line)
                if (only_resized_triangle(j_component)%line(i_component)) then
                    call this%short_interactions%cells(i_box)%neighbour_cells(j_component)%line(i_component)%&
                        cells%reset()
                else
                    call cells_destroy(this%short_interactions%cells(i_box)%neighbour_cells(j_component)%&
                        line(i_component)%cells)
                    allocate(this%short_interactions%cells(i_box)%neighbour_cells(j_component)%&
                        line(i_component)%cells, source=neighbour_cells(j_component)%&
                        line(i_component)%cells)
                end if
            end do
        end do
        call cells_destroy(neighbour_cells)
        call this%short_interactions%visitable_cells_memento%restore(this%short_interactions%cells(i_box)%&
            visitable_cells, this%short_interactions%cells(i_box)%neighbour_cells, only_resized_triangle, &
            visitable_cells)
        call cells_destroy(visitable_cells)
    end subroutine Abstract_restore_cells

end module classes_box_volume_change
