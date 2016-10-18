module classes_volume_change_method

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive
use types_logical_line, only: Concrete_Logical_Line
use procedures_logical_factory, only: logical_create => create
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
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
use classes_changed_box_size_ratio, only: Abstract_Changed_Box_Size_Ratio
use procedures_triangle_observables, only: operator(-)
use types_observables_energies, only: Concrete_Observables_Energies
use procedures_triangle_observables, only: triangle_observables_sum
use procedures_observables_energies_factory, only: observables_energies_create => create
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use classes_exploring_algorithm, only: Abstract_Exploring_Algorithm

implicit none

private

    type, extends(Abstract_Exploring_Algorithm), abstract, public :: Abstract_Volume_Change_Method
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Component_Wrapper), pointer :: components(:, :) => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        class(Abstract_Dipolar_Interactions_Facade), pointer :: dipolar_interactions_facades(:) => &
            null()
        class(Abstract_Changed_Box_Size_Ratio), pointer :: changed_boxes_size_ratio(:) => null()
        integer :: num_changes = 0
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
        dipolar_interactions_facades, changed_boxes_size_ratio, num_changes)
        class(Abstract_Volume_Change_Method), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:, :)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        class(Abstract_Dipolar_Interactions_Facade), target, intent(in) :: &
            dipolar_interactions_facades(:)
        class(Abstract_Changed_Box_Size_Ratio), target, intent(in) :: changed_boxes_size_ratio(:)
        integer, intent(in) :: num_changes

        this%environment => environment
        this%components => components
        this%short_interactions => short_interactions
        this%dipolar_interactions_facades => dipolar_interactions_facades
        this%changed_boxes_size_ratio => changed_boxes_size_ratio
        call check_positive("Abstract_Volume_Change_Method: construct", "num_changes", num_changes)
        this%num_changes = num_changes
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Volume_Change_Method), intent(inout) :: this

        this%changed_boxes_size_ratio => null()
        this%dipolar_interactions_facades => null()
        this%short_interactions => null()
        this%components => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_try(this, observables)
        class(Abstract_Volume_Change_Method), intent(in) :: this
        type(Exploring_Observables_Wrapper), intent(inout) :: observables !too much?

        logical :: overlap
        integer :: i_box
        real(DP) :: beta_pressure_excess_sum
        real(DP) :: contacts, delta_energy, delta_volume
        real(DP), dimension(num_dimensions) :: new_box_size, box_size, box_size_ratio
        type(Neighbour_Cells_Line), allocatable :: neighbour_cells(:)
        type(Concrete_Logical_Line), allocatable :: only_resized_triangle(:)
        class(Abstract_Visitable_Cells), allocatable :: visitable_cells(:, :)
        type(Dipolar_Interactions_Static_Wrapper) :: dipolar_interactions_static
        logical :: reset_real_pair
        integer :: i_change

        do i_box = 1, size(this%environment%periodic_boxes)
            call short_interactions_visit_cells(overlap, contacts, this%components(:, i_box), this%&
                short_interactions%cells(i_box)%visitable_cells)
            if (overlap) then
                call error_exit("Abstract_Volume_Change_Method: try: "//&
                    "short_interactions_visit_cells: overlap")
            end if

            box_size = this%environment%periodic_boxes(i_box)%get_size()
            beta_pressure_excess_sum = 0._DP
            do i_change = 1, this%num_changes
                box_size_ratio = this%changed_boxes_size_ratio(i_box)%get()
                new_box_size = box_size * box_size_ratio
                call this%dipolar_interactions_facades(i_box)%save(dipolar_interactions_static, &
                    reset_real_pair, new_box_size)

                call this%environment%periodic_boxes(i_box)%set(box_size * box_size_ratio)
                call this%rescale_positions(i_box, box_size_ratio)
                call logical_create(only_resized_triangle, size(this%components, 1))
                call this%save_cells(neighbour_cells, only_resized_triangle, visitable_cells, i_box)
                call short_interactions_reset(this%short_interactions%cells(i_box)%neighbour_cells,&
                    only_resized_triangle, this%short_interactions%cells(i_box)%visitable_cells)
                call this%dipolar_interactions_facades(i_box)%reset()
                call this%set_delta_energy(overlap, delta_energy, i_box, box_size_ratio, &
                    observables%energies(i_box))
                if (overlap) call error_exit("Abstract_Volume_Change_Method: try: "//&
                    "set_delta_energy: overlap")
                delta_volume = product(this%environment%accessible_domains(i_box)%get_size()) * &
                    (product(box_size_ratio) - 1._DP)
                beta_pressure_excess_sum = beta_pressure_excess_sum - delta_energy / delta_volume /&
                    this%environment%temperature%get()

                call this%environment%periodic_boxes(i_box)%set(box_size)
                call this%rescale_positions(i_box, 1._DP / box_size_ratio)
                call this%restore_cells(neighbour_cells, only_resized_triangle, visitable_cells, &
                    i_box)
                call this%dipolar_interactions_facades(i_box)%restore(dipolar_interactions_static, &
                    reset_real_pair)
            end do
            observables%beta_pressures_excess(i_box) = this%short_interactions%&
                beta_pressures_excess(i_box)%get(contacts) + &
                beta_pressure_excess_sum / real(this%num_changes, DP)
        end do
    end subroutine Abstract_try


    subroutine Abstract_set_delta_energy(this, overlap, delta_energy, i_box, box_size_ratio, &
        energies)
        class(Abstract_Volume_Change_Method), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_box
        real(DP), intent(in) :: box_size_ratio(:)
        type(Concrete_Observables_Energies), intent(in) :: energies

        type(Concrete_Observables_Energies) :: deltas, new_energies

        call observables_energies_create(deltas, size(this%components, 1))
        call observables_energies_create(new_energies, size(this%components, 1))
        call short_interactions_visit(overlap, new_energies%walls_energies, this%&
            components(:, i_box), this%short_interactions%walls_visitors(i_box), this%&
            short_interactions%wall_pairs)
        if (overlap) return
        deltas%walls_energies = new_energies%walls_energies - energies%walls_energies
        call short_interactions_visit_cells(overlap, new_energies%short_energies, &
            this%components(:, i_box), this%short_interactions%cells(i_box)%visitable_cells)
        if (overlap) return
        deltas%short_energies = new_energies%short_energies - energies%short_energies
        call dipolar_interactions_visit(new_energies%field_energies, this%environment%&
            external_fields(i_box), this%components(:, i_box))
        deltas%field_energies = new_energies%field_energies - energies%field_energies
        call this%dipolar_interactions_facades(i_box)%visit(new_energies%dipolar_energies, &
            new_energies%dipolar_shared_energy, product(box_size_ratio), energies%dipolar_energies,&
            energies%dipolar_shared_energy)
        deltas%dipolar_energies = new_energies%dipolar_energies - energies%dipolar_energies
        deltas%dipolar_shared_energy = new_energies%dipolar_shared_energy - energies%&
            dipolar_shared_energy

        delta_energy = sum(deltas%walls_energies + deltas%field_energies) + &
            triangle_observables_sum(deltas%short_energies) + &
            triangle_observables_sum(deltas%dipolar_energies) + deltas%dipolar_shared_energy
    end subroutine Abstract_set_delta_energy

    subroutine Abstract_rescale_positions(this, i_box, box_size_ratio)
        class(Abstract_Volume_Change_Method), intent(in) :: this
        integer, intent(in) :: i_box
        real(DP), intent(in) :: box_size_ratio(:)

        integer :: i_component

        do i_component = 1, size(this%components, 1)
            call this%components(i_component, i_box)%positions%rescale_all(box_size_ratio)
        end do
    end subroutine Abstract_rescale_positions

    !> @note
    !> cf. [[classes_box_volume_change:Abstract_save_cells]]
    subroutine Abstract_save_cells(this, neighbour_cells, only_resized_triangle, visitable_cells, &
        i_box)
        class(Abstract_Volume_Change_Method), intent(in) :: this
        type(Neighbour_Cells_Line), allocatable, intent(out) :: neighbour_cells(:)
        type(Concrete_Logical_Line), intent(inout) ::only_resized_triangle(:)
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: visitable_cells(:, :)
        integer, intent(in) :: i_box

        integer :: i_component, j_component

        call cells_allocate_triangle(neighbour_cells, size(this%components, 1))
        do j_component = 1, size(neighbour_cells)
            do i_component = 1, size(neighbour_cells(j_component)%line)
                only_resized_triangle(j_component)%line(i_component) = this%short_interactions%&
                    cells(i_box)%neighbour_cells(j_component)%line(i_component)%cells%resize_only()
                if (.not. only_resized_triangle(j_component)%line(i_component)) then
                    allocate(neighbour_cells(j_component)%line(i_component)%cells, source=this%&
                        short_interactions%cells(i_box)%neighbour_cells(j_component)%&
                        line(i_component)%cells)
                end if
            end do
        end do
        call this%short_interactions%visitable_cells_memento%save(visitable_cells, this%&
            short_interactions%cells(i_box)%visitable_cells)
    end subroutine Abstract_save_cells

    !> @note
    !> cf. [[classes_box_volume_change:Abstract_restore_cells]]
    subroutine Abstract_restore_cells(this, neighbour_cells, only_resized_triangle, &
        visitable_cells, i_box)
        class(Abstract_Volume_Change_Method), intent(in) :: this
        type(Neighbour_Cells_Line), allocatable, intent(inout) :: neighbour_cells(:)
        type(Concrete_Logical_Line), intent(in) ::only_resized_triangle(:)
        class(Abstract_Visitable_Cells), allocatable, intent(inout) :: visitable_cells(:, :)
        integer, intent(in) :: i_box

        integer :: i_component, j_component

        do j_component = 1, size(this%short_interactions%cells(i_box)%neighbour_cells)
            do i_component = 1, size(this%short_interactions%cells(i_box)%&
                neighbour_cells(j_component)%line)
                if (only_resized_triangle(j_component)%line(i_component)) then
                    call this%short_interactions%cells(i_box)%neighbour_cells(j_component)%&
                        line(i_component)%cells%reset()
                else
                    call cells_destroy(this%short_interactions%cells(i_box)%&
                        neighbour_cells(j_component)%line(i_component)%cells)
                    allocate(this%short_interactions%cells(i_box)%neighbour_cells(j_component)%&
                        line(i_component)%cells, source=neighbour_cells(j_component)%&
                        line(i_component)%cells)
                end if
            end do
        end do
        call cells_destroy(neighbour_cells)
        call this%short_interactions%visitable_cells_memento%restore(this%short_interactions%&
            cells(i_box)%visitable_cells, this%short_interactions%cells(i_box)%&
            neighbour_cells, only_resized_triangle, visitable_cells)
        call cells_destroy(visitable_cells)
    end subroutine Abstract_restore_cells

!end implementation Abstract_Volume_Change_Method

!implementation Null_Volume_Change_Method

    subroutine Null_construct(this, environment, components, short_interactions, &
        dipolar_interactions_facades, changed_boxes_size_ratio, num_changes)
        class(Null_Volume_Change_Method), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:, :)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        class(Abstract_Dipolar_Interactions_Facade), target, intent(in) :: &
            dipolar_interactions_facades(:)
        class(Abstract_Changed_Box_Size_Ratio), target, intent(in) :: changed_boxes_size_ratio(:)
        integer, intent(in) :: num_changes
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
