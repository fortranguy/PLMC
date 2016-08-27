module classes_volume_change_method

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_errors, only: error_exit
use types_logical_line, only: Concrete_Logical_Line
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_neighbour_cells_wrapper, only: Neighbour_Cells_Line
use classes_visitable_cells, only: Abstract_Visitable_Cells
use procedures_cells_factory, only: cells_destroy => destroy
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use classes_changed_box_size_ratio, only: Abstract_Changed_Box_Size_Ratio
use procedures_changed_box_size_ratio_factory, only: changed_box_size_ratio_destroy => destroy
use procedures_triangle_observables, only: operator(-)
use types_observables_energies, only: Concrete_Energies
use procedures_triangle_observables, only: triangle_observables_sum
use procedures_observables_energies_factory, only: observables_energies_create => create
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use classes_exploring_algorithm, only: Abstract_Exploring_Algorithm
use procedures_plmc_reset, only: box_size_change_reset_cells
use procedures_plmc_visit, only: visit_walls, visit_short, visit_field

implicit none

private

    type, extends(Abstract_Exploring_Algorithm), abstract, public :: Abstract_Volume_Change_Method
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Component_Wrapper), pointer :: components(:) => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        class(Abstract_Changed_Box_Size_Ratio), allocatable :: changed_box_size_ratio
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: try => Abstract_try
        procedure, private :: set_delta_energy => Abstract_set_delta_energy
        procedure, private :: rescale_positions => Abstract_rescale_positions
        procedure, private :: save_cells => Abstract_save_cells
        procedure, private :: restore_cells => Abstract_restore_cells
    end type Abstract_Volume_Change_Method

    type, extends(Abstract_Volume_Change_Method), public :: Concrete_Volume_Change_Method

    end type Concrete_Volume_Change_Method

    type, extends(Abstract_Volume_Change_Method), public :: Null_Volume_Change_Method
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: try => Null_try
    end type Null_Volume_Change_Method

contains

!implementation Abstract_Volume_Change_Method

    subroutine Abstract_construct(this, environment, components, short_interactions, &
        changed_box_size_ratio)
        class(Abstract_Volume_Change_Method), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        class(Abstract_Changed_Box_Size_Ratio), intent(in) :: changed_box_size_ratio

        this%environment => environment
        this%components => components
        this%short_interactions => short_interactions
        allocate(this%changed_box_size_ratio, source=changed_box_size_ratio)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Volume_Change_Method), intent(inout) :: this

        call changed_box_size_ratio_destroy(this%changed_box_size_ratio)
        this%short_interactions => null()
        this%components => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_try(this, observables)
        class(Abstract_Volume_Change_Method), intent(in) :: this
        type(Exploring_Observables_Wrapper), intent(inout) :: observables !too much?

        logical :: overlap
        real(DP) :: contacts, delta_energy, d_energy_d_volume
        real(DP), dimension(num_dimensions) :: box_size, box_size_ratio
        type(Neighbour_Cells_Line), allocatable :: neighbour_cells(:)
        type(Concrete_Logical_Line), allocatable :: only_resized_triangle(:)
        class(Abstract_Visitable_Cells), allocatable :: visitable_cells(:, :)

        call visit_short(overlap, contacts, this%components, this%short_interactions)
        if (overlap) call error_exit("Abstract_Volume_Change_Method: try: visit_short: overlap")
        call this%save_cells(neighbour_cells, visitable_cells)
        box_size_ratio = this%changed_box_size_ratio%get()
        box_size = this%environment%periodic_box%get_size()
        call this%environment%periodic_box%set(box_size * box_size_ratio)
        call this%rescale_positions(box_size_ratio)
        call box_size_change_reset_cells(this%short_interactions%neighbour_cells, &
            only_resized_triangle, this%short_interactions%visitable_cells)
        call this%set_delta_energy(overlap, delta_energy, observables%energies)
        if (overlap) call error_exit("Abstract_Volume_Change_Method: try: set_delta_energy: "//&
            "overlap")

        d_energy_d_volume = delta_energy / product(this%environment%accessible_domain%get_size()) /&
            (product(box_size_ratio) - 1._DP)
        observables%beta_pressure_excess = this%short_interactions%beta_pressure_excess%&
            get(contacts) - d_energy_d_volume / this%environment%temperature%get()

        call this%environment%periodic_box%set(box_size)
        call this%rescale_positions(1._DP / box_size_ratio)
        call this%restore_cells(neighbour_cells, only_resized_triangle, visitable_cells)
    end subroutine Abstract_try

    subroutine Abstract_set_delta_energy(this, overlap, delta_energy, energies)
        class(Abstract_Volume_Change_Method), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energy
        type(Concrete_Energies), intent(in) :: energies

        type(Concrete_Energies) :: deltas, new_energies

        call observables_energies_create(deltas, size(this%components))
        call observables_energies_create(new_energies, size(this%components))
        call visit_walls(overlap, new_energies%walls_energies, this%components, this%&
            short_interactions)
        if (overlap) return
        deltas%walls_energies = new_energies%walls_energies - energies%walls_energies
        call visit_short(overlap, new_energies%short_energies, this%components, this%&
            short_interactions)
        if (overlap) return
        deltas%short_energies = new_energies%short_energies - energies%short_energies
        call visit_field(new_energies%field_energies, this%environment%external_field, this%&
            components)
        deltas%field_energies = new_energies%field_energies - energies%field_energies
        !dipolar interactions

        delta_energy = sum(deltas%walls_energies + deltas%field_energies) + &
            triangle_observables_sum(deltas%short_energies)
    end subroutine Abstract_set_delta_energy

    subroutine Abstract_rescale_positions(this, box_size_ratio)
        class(Abstract_Volume_Change_Method), intent(in) :: this
        real(DP), intent(in) :: box_size_ratio(:)

        integer :: i_component

        do i_component = 1, size(this%components)
            call this%components(i_component)%positions%rescale_all(box_size_ratio)
        end do
    end subroutine Abstract_rescale_positions

    subroutine Abstract_save_cells(this, neighbour_cells, visitable_cells)
        class(Abstract_Volume_Change_Method), intent(in) :: this
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
        call this%short_interactions%visitable_cells_memento%save(visitable_cells, this%&
            short_interactions%visitable_cells)
    end subroutine Abstract_save_cells

    subroutine Abstract_restore_cells(this, neighbour_cells, only_resized_triangle, visitable_cells)
        class(Abstract_Volume_Change_Method), intent(in) :: this
        type(Neighbour_Cells_Line), intent(in) :: neighbour_cells(:)
        type(Concrete_Logical_Line), intent(in) ::only_resized_triangle(:)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells(:, :)

        integer :: i_component, j_component

        call cells_destroy(this%short_interactions%neighbour_cells)
        allocate(this%short_interactions%neighbour_cells(size(neighbour_cells)))
        do j_component = 1, size(neighbour_cells)
            allocate(this%short_interactions%neighbour_cells(j_component)%line(j_component))
            do i_component = 1, size(neighbour_cells(j_component)%line)
                allocate(this%short_interactions%neighbour_cells(j_component)%line(i_component)%&
                    cells, source=neighbour_cells(j_component)%line(i_component)%cells)
            end do
        end do
        call this%short_interactions%visitable_cells_memento%restore(this%short_interactions%&
            visitable_cells, this%short_interactions%neighbour_cells, only_resized_triangle, &
            visitable_cells)
    end subroutine Abstract_restore_cells

!end implementation Abstract_Volume_Change_Method

!implementation Null_Volume_Change_Method

    subroutine Null_construct(this, environment, components, short_interactions, &
        changed_box_size_ratio)
        class(Null_Volume_Change_Method), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        class(Abstract_Changed_Box_Size_Ratio), intent(in) :: changed_box_size_ratio
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Volume_Change_Method), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_try(this, observables)
        class(Null_Volume_Change_Method), intent(in) :: this
        type(Exploring_Observables_Wrapper), intent(inout) :: observables
    end subroutine Null_try

!end implementation Null_Volume_Change_Method

end module classes_volume_change_method
