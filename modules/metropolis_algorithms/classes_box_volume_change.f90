module classes_box_volume_change

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_neighbour_cells_wrapper, only: Neighbour_Cells_Line
use classes_visitable_cells, only: Abstract_Visitable_Cells
use procedures_cells_factory, only: cells_destroy => destroy
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_visit_condition, only: visit_condition_in_range => in_range, &
    visit_condition_lower => lower, visit_condition_unconditional => unconditional
use procedures_dipoles_field_interaction, only: dipoles_field_visit_component => visit_component
use classes_changed_box_size, only: Abstract_Changed_Box_Size
use types_reals_line, only: Reals_Line
use procedures_triangle_observables, only: operator(-)
use types_observables_energies, only: Concrete_Single_Energies, Concrete_Energies
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use procedures_triangle_observables, only: triangle_observables_sum
use procedures_observables_energies_factory, only: observables_energies_create => create
use procedures_generating_observables_factory, only:generating_observables_create => create
use classes_metropolis_algorithm, only: Abstract_Metropolis_Algorithm
use procedures_plmc_reset, only: plmc_reset_cells, plmc_reset_cells

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
        procedure, private :: test_metropolis => Abstract_test_metropolis
        procedure, private :: acceptation_probability => Abstract_acceptation_probability
        procedure, private :: rescale_positions => Abstract_rescale_positions
        procedure, private :: save_cells => Abstract_save_cells
        procedure, private :: restore_cells => Abstract_restore_cells
        procedure, private :: visit_walls => Abstrac_visit_walls
        procedure, private :: visit_fields => Abstract_visit_fields
        procedure, private :: visit_short => Abstract_visit_short
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
        call observables_energies_create(deltas, size(this%mixture%components))
        call this%test_metropolis(success, deltas, observables%energies)
    end subroutine Abstract_try

    subroutine Abstract_test_metropolis(this, success, deltas, energies)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Energies), intent(inout) :: deltas
        type(Concrete_Energies), intent(in) :: energies

        real(DP) :: delta_energy
        type(Concrete_Energies) :: new_energies
        real(DP), dimension(num_dimensions) :: box_size, box_size_ratio
        type(Neighbour_Cells_Line), allocatable :: neighbour_cells(:)
        class(Abstract_Visitable_Cells), allocatable :: visitable_cells(:, :)
        logical :: overlap
        real(DP) :: acceptation_probability, rand

        success = .false.
        box_size_ratio = this%changed_box_size%get_ratio()
        box_size = this%environment%periodic_box%get_size()
        call this%environment%periodic_box%set(box_size * box_size_ratio)
        call this%rescale_positions(box_size_ratio)
        call this%save_cells(neighbour_cells, visitable_cells)
        call plmc_reset_cells(this%short_interactions%neighbour_cells)
        call plmc_reset_cells(this%short_interactions%visitable_cells)

        call generating_observables_create(new_energies%walls_energies, size(deltas%walls_energies))
        call this%visit_walls(overlap, new_energies%walls_energies)
        if (overlap) return
        deltas%walls_energies = new_energies%walls_energies - energies%walls_energies
        call generating_observables_create(new_energies%field_energies, size(deltas%field_energies))
        call this%visit_fields(new_energies%field_energies)
        deltas%field_energies = new_energies%field_energies - energies%field_energies
        call generating_observables_create(new_energies%short_energies, size(deltas%short_energies))
        call this%visit_short(overlap, new_energies%short_energies)
        if (overlap) return
        deltas%short_energies = new_energies%short_energies - energies%short_energies
        !dipolar interactions

        delta_energy = sum(deltas%walls_energies + deltas%field_energies) + &
            triangle_observables_sum(deltas%short_energies)
        acceptation_probability = this%acceptation_probability(box_size_ratio, delta_energy)
        if (acceptation_probability < 1._DP) then
            call random_number(rand)
            if (rand < acceptation_probability) success = .true.
        else
            success = .true.
        end if

        if (.not.success) then
            call this%environment%periodic_box%set(box_size)
            call this%rescale_positions(1._DP / box_size_ratio)
            call this%restore_cells(neighbour_cells, visitable_cells)
        end if
    end subroutine Abstract_test_metropolis

    !> \[
    !>      P[V \to V^\prime] = e^{-\beta p V \left( \frac{V^\prime}{V} - 1 \right)}
    !>          \left( \frac{V^\prime}{V} \right)^{N+1}
    !>          e^{-\beta [U(\vec{s}^N, V^\prime) - U(\vec{s}^N, V)]}
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
            exp(-delta_energy/this%environment%temperature%get())
        if (probability > 1._DP) probability = 1._DP
    end function Abstract_acceptation_probability

        subroutine Abstract_rescale_positions(this, box_size_ratio)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        real(DP), intent(in) :: box_size_ratio(:)

        integer :: i_component

        do i_component = 1, size(this%mixture%components)
            call this%mixture%components(i_component)%positions%rescale_all(box_size_ratio)
        end do
    end subroutine Abstract_rescale_positions

    subroutine Abstract_save_cells(this, neighbour_cells, visitable_cells)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        type(Neighbour_Cells_Line), allocatable, intent(out) :: neighbour_cells(:)
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: visitable_cells(:, :)

        integer :: i_component, j_component

        allocate(neighbour_cells(size(this%short_interactions%neighbour_cells)))
        do j_component = 1, size(neighbour_cells)
            allocate(neighbour_cells(j_component)%line(j_component))
            do i_component = 1, size(neighbour_cells(j_component)%line)
                allocate(neighbour_cells(j_component)%line(i_component)%cells, source=this%&
                    short_interactions%neighbour_cells(j_component)%line(i_component)%cells)
            end do
        end do
        allocate(visitable_cells, source=this%short_interactions%visitable_cells)
    end subroutine Abstract_save_cells

    subroutine Abstract_restore_cells(this, neighbour_cells, visitable_cells)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        type(Neighbour_Cells_Line), intent(in) :: neighbour_cells(:)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells(:, :)

        integer :: i_component, j_component

        call cells_destroy(this%short_interactions%visitable_cells)
        call cells_destroy(this%short_interactions%neighbour_cells)

        do j_component = 1, size(neighbour_cells)
            do i_component = 1, size(neighbour_cells(j_component)%line)
                allocate(this%short_interactions%neighbour_cells(j_component)%line(i_component)%&
                    cells, source=neighbour_cells(j_component)%line(i_component)%cells)
            end do
        end do
        allocate(this%short_interactions%visitable_cells, source=visitable_cells)
    end subroutine Abstract_restore_cells

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
                components(i_component)%dipole_moments)
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

!end implementation Abstract_Box_Volume_Change

!implementation Null_Box_Volume_Change

    subroutine Null_construct(this, environment, mixture, short_interactions, changed_box_size)
        class(Null_Box_Volume_Change), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
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
