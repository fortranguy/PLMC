module classes_box_volume_change

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use types_logical_line, only: Concrete_Logical_Line
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_neighbour_cells_wrapper, only: Neighbour_Cells_Line
use classes_visitable_cells, only: Abstract_Visitable_Cells
use procedures_cells_memento, only: cells_save => save, cells_restore => restore
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_short_interactions_resetter, only: short_interactions_reset => reset
use procedures_short_interactions_visitor, only: short_interactions_visit => visit, &
    short_interactions_visit_cells => visit_cells
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use procedures_dipolar_interactions_visitor, only: dipolar_interactions_visit => visit
use classes_dipolar_interactions_facade, only: Abstract_Dipolar_Interactions_Facade
use classes_changed_box_size, only: Abstract_Changed_Box_Size
use procedures_triangle_observables, only: operator(-)
use types_observables_energies, only: Concrete_Energies
use procedures_observables_energies_factory, only: observables_energies_set => set
use procedures_triangle_observables, only: triangle_observables_sum
use procedures_observables_energies_factory, only: observables_energies_create => create
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use classes_generating_algorithm, only: Abstract_Generating_Algorithm
use procedures_metropolis_algorithm, only: metropolis_algorithm

implicit none

private

    type, extends(Abstract_Generating_Algorithm), abstract, public :: Abstract_Box_Volume_Change
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        class(Abstract_Dipolar_Interactions_Facade), pointer :: dipolar_interactions_facade => &
            null()
        class(Abstract_Changed_Box_Size), pointer :: changed_box_size => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: get_num_choices => Abstract_get_num_choices
        procedure :: try => Abstract_try
        procedure, private :: metropolis_algorithm => Abstract_metropolis_algorithm
        procedure, private :: acceptation_probability => Abstract_acceptation_probability
        procedure, private :: rescale_positions => Abstract_rescale_positions
    end type Abstract_Box_Volume_Change

    type, extends(Abstract_Box_Volume_Change), public :: Concrete_Box_Volume_Change

    end type Concrete_Box_Volume_Change

    type, extends(Abstract_Box_Volume_Change), public :: Null_Box_Volume_Change
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: get_num_choices => Null_get_num_choices
        procedure :: try => Null_try
    end type Null_Box_Volume_Change

contains

!implementation Abstract_Box_Volume_Change

    subroutine Abstract_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_facade, changed_box_size)
        class(Abstract_Box_Volume_Change), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        class(Abstract_Dipolar_Interactions_Facade), target, intent(in) :: &
            dipolar_interactions_facade
        class(Abstract_Changed_Box_Size), target, intent(in) :: changed_box_size

        this%environment => environment
        this%mixture => mixture
        this%short_interactions => short_interactions
        this%dipolar_interactions_facade => dipolar_interactions_facade
        this%changed_box_size => changed_box_size
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Box_Volume_Change), intent(inout) :: this

        this%changed_box_size => null()
        this%dipolar_interactions_facade => null()
        this%short_interactions => null()
        this%mixture => null()
        this%environment => null()
    end subroutine Abstract_destroy

    pure integer function Abstract_get_num_choices(this) result(num_choices)
        class(Abstract_Box_Volume_Change), intent(in) :: this

        num_choices = this%changed_box_size%get_num()
    end function Abstract_get_num_choices

    subroutine Abstract_try(this, observables)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Energies) :: new_energies
        real(DP), dimension(num_dimensions) :: new_box_size, box_size, box_size_ratio
        type(Neighbour_Cells_Line), allocatable :: neighbour_cells(:)
        type(Concrete_Logical_Line), allocatable :: only_resized_triangle(:)
        class(Abstract_Visitable_Cells), allocatable :: visitable_cells(:, :)
        type(Dipolar_Interactions_Static_Wrapper) :: dipolar_interactions_static

        observables%volume_change_counter%num_hits = observables%volume_change_counter%num_hits + 1
        box_size = this%environment%periodic_box%get_size()
        box_size_ratio = this%changed_box_size%get_ratio()
        new_box_size = box_size * box_size_ratio
        call cells_save(this%short_interactions%visitable_cells_memento, neighbour_cells, &
            visitable_cells, this%short_interactions%neighbour_cells, this%short_interactions%&
            visitable_cells)
        call this%dipolar_interactions_facade%save(dipolar_interactions_static, &
            product(new_box_size))

        call this%environment%periodic_box%set(new_box_size)
        call this%rescale_positions(box_size_ratio)
        call short_interactions_reset(this%short_interactions%neighbour_cells, &
            only_resized_triangle, this%short_interactions%visitable_cells)
        call this%dipolar_interactions_facade%reset()
        call observables_energies_create(new_energies, size(this%mixture%components))
        call this%metropolis_algorithm(success, new_energies, box_size_ratio, observables%energies)

        if (success) then
            observables%accessible_domain_size = this%environment%accessible_domain%get_size()
            call observables_energies_set(observables%energies, new_energies)
            observables%volume_change_counter%num_successes = observables%volume_change_counter%&
                num_successes + 1
        else
            call this%environment%periodic_box%set(box_size)
            call this%rescale_positions(1._DP / box_size_ratio)
            call cells_restore(this%short_interactions%visitable_cells_memento, this%&
                short_interactions%neighbour_cells, this%short_interactions%visitable_cells, &
                neighbour_cells, only_resized_triangle, visitable_cells)
            call this%dipolar_interactions_facade%restore(dipolar_interactions_static)
        end if
    end subroutine Abstract_try

    subroutine Abstract_metropolis_algorithm(this, success, new_energies, box_size_ratio, energies)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Energies), intent(inout) :: new_energies
        real(DP), intent(in) :: box_size_ratio(:)
        type(Concrete_Energies), intent(in) :: energies

        real(DP) :: delta_energy
        type(Concrete_Energies) :: deltas

        logical :: overlap

        success = .false.
        call observables_energies_create(deltas, size(this%mixture%components))
        call short_interactions_visit(overlap, new_energies%walls_energies, this%mixture%&
            components, this%short_interactions)
        if (overlap) return
        deltas%walls_energies = new_energies%walls_energies - energies%walls_energies
        call short_interactions_visit_cells(overlap, new_energies%short_energies, this%mixture%&
            components, this%short_interactions)
        deltas%short_energies = new_energies%short_energies - energies%short_energies
        if (overlap) return
        call dipolar_interactions_visit(new_energies%field_energies, this%environment%&
            external_field, this%mixture%components)
        deltas%field_energies = new_energies%field_energies - energies%field_energies
        call this%dipolar_interactions_facade%visit(new_energies%dipolar_energies, new_energies%&
            dipolar_shared_energy, product(box_size_ratio), energies%dipolar_energies, energies%&
            dipolar_shared_energy)
        deltas%dipolar_energies = new_energies%dipolar_energies - energies%dipolar_energies
        deltas%dipolar_shared_energy = new_energies%dipolar_shared_energy - energies%&
            dipolar_shared_energy

        delta_energy = sum(deltas%walls_energies + deltas%field_energies) + &
            triangle_observables_sum(deltas%short_energies) + &
            triangle_observables_sum(deltas%dipolar_energies) + deltas%dipolar_shared_energy
        success = metropolis_algorithm(this%acceptation_probability(box_size_ratio, delta_energy))
    end subroutine Abstract_metropolis_algorithm

    !> \[
    !>      P[V \to V^\prime] = min \left( 1, e^{-\beta p V \left( \frac{V^\prime}{V} - 1 \right)}
    !>          \left( \frac{V^\prime}{V} \right)^{N+1}
    !>          e^{-\beta [U(\vec{s}^N, V^\prime) - U(\vec{s}^N, V)]} \right)
    !> \]
    pure real(DP) function Abstract_acceptation_probability(this, box_size_ratio, delta_energy) &
        result(probability)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        real(DP), intent(in) :: box_size_ratio(:)
        real(DP), intent(in) :: delta_energy

        real(DP) :: volume_ratio
        integer :: i_component, num_particles

        num_particles = 0._DP
        do i_component = 1, size(this%mixture%components)
            num_particles = num_particles + this%mixture%components(i_component)%number%get()
        end do
        volume_ratio = product(box_size_ratio)

        probability = exp(-this%environment%beta_pressure%get() * &
            product(this%environment%accessible_domain%get_size()) * (volume_ratio - 1._DP)) * &
            volume_ratio**(num_particles + 1) * &
            exp(-delta_energy / this%environment%temperature%get())
        probability = min(1._DP, probability)
    end function Abstract_acceptation_probability

    subroutine Abstract_rescale_positions(this, box_size_ratio)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        real(DP), intent(in) :: box_size_ratio(:)

        integer :: i_component

        do i_component = 1, size(this%mixture%components)
            call this%mixture%components(i_component)%positions%rescale_all(box_size_ratio)
        end do
    end subroutine Abstract_rescale_positions

!end implementation Abstract_Box_Volume_Change

!implementation Null_Box_Volume_Change

    subroutine Null_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_facade, changed_box_size)
        class(Null_Box_Volume_Change), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        class(Abstract_Dipolar_Interactions_Facade), target, intent(in) :: &
            dipolar_interactions_facade
        class(Abstract_Changed_Box_Size), target, intent(in) :: changed_box_size
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Box_Volume_Change), intent(inout) :: this
    end subroutine Null_destroy

    pure integer function Null_get_num_choices(this) result(num_choices)
        class(Null_Box_Volume_Change), intent(in) :: this
        num_choices = 0
    end function Null_get_num_choices

    subroutine Null_try(this, observables)
        class(Null_Box_Volume_Change), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables
    end subroutine Null_try

!end implementation Null_Box_Volume_Change

end module classes_box_volume_change
