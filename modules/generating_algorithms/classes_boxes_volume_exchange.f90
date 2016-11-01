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
use types_mixture_wrapper, only: Mixture_Wrapper
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
use classes_exchanged_boxes_size, only: Abstract_Exchanged_Boxes_Size
use procedures_triangle_observables, only: triangle_observables_diff, triangle_observables_sum
use types_observables_energies, only: Concrete_Observables_Energies
use procedures_observables_energies_factory, only: observables_energies_create => create, &
    observables_energies_set => set
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use classes_generating_algorithm, only: Abstract_Generating_Algorithm
use procedures_selectors_resetters, only: selectors_reset => reset

implicit none

private

    type, extends(Abstract_Generating_Algorithm), public :: Boxes_Volume_Exchange
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        class(Abstract_Dipolar_Interactions_Facade), pointer :: dipolar_interactions_facades(:) => &
            null()
        class(Abstract_Exchanged_Boxes_Size), pointer :: exchanged_boxes_size => null()
        logical, allocatable :: have_positions(:, :)
        class(Abstract_Hetero_Couples), allocatable :: couples
        class(Abstract_Tower_Sampler), allocatable :: selector
    contains
        procedure :: construct => Concrete_construct
        procedure :: destroy => Concrete_destroy
        procedure :: reset_selectors => Concrete_reset_selectors
        procedure :: get_num_choices => Concrete_get_num_choices
        procedure :: try => Concrete_try
    end type Boxes_Volume_Exchange

contains

    subroutine Concrete_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_facades, exchanged_boxes_size, have_positions, couples, selector)
        class(Boxes_Volume_Exchange), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        class(Abstract_Dipolar_Interactions_Facade), target, intent(in) :: &
            dipolar_interactions_facades(:)
        class(Abstract_Exchanged_Boxes_Size), target, intent(in) :: exchanged_boxes_size
        logical, intent(in) :: have_positions(:, :)
        class(Abstract_Hetero_Couples), intent(in) :: couples
        class(Abstract_Tower_Sampler), intent(in) :: selector

        this%environment => environment
        this%mixture => mixture
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
        this%mixture => null()
        this%environment => null()
    end subroutine Concrete_destroy

    subroutine Concrete_reset_selectors(this)
        class(Boxes_Volume_Exchange), intent(inout) :: this

        call selectors_reset(this%selector, this%couples, this%exchanged_boxes_size, this%mixture%&
            components, this%have_positions)
    end subroutine Concrete_reset_selectors

    pure integer function Concrete_get_num_choices(this) result(num_choices)
        class(Boxes_Volume_Exchange), intent(in) :: this

        num_choices = this%selector%get_num_choices()
    end function Concrete_get_num_choices

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

        !observables%volumes_change_counter(j_box)%line(i_box)%num_hits = observables%&
        !    volumes_change_counter(j_box)%line(i_box)%num_hits + 1

        do i_partner = 1, size(new_boxes_size, 2)
            boxes_size(:, i_partner) = this%environment%periodic_boxes(ij_boxes(i_partner))%&
                get_size()
        end do
        boxes_size_ratio = this%exchanged_boxes_size%get_ratios(boxes_size(:, 1) / boxes_size(:, 2))
        new_boxes_size = boxes_size * boxes_size_ratio

        do i_partner = 1, size(dipolar_interactions_static)
            call this%dipolar_interactions_facades(ij_boxes(i_partner))%&
                save(dipolar_interactions_static(i_partner), reset_real_pair(i_partner), &
                    new_boxes_size(:, i_partner))
        end do

        do i_partner = 1, size(new_boxes_size, 2)
            call this%environment%periodic_boxes(ij_boxes(i_partner))%&
                set(new_boxes_size(:, i_partner))
            call mixture_rescale_positions(this%mixture%components(:, ij_boxes(i_partner)), &
                boxes_size_ratio(:, i_partner))
            call logical_create(only_resized_triangle(i_partner)%triangle, &
                size(this%mixture%components, 1))
            call cells_memento_save(cells(i_partner), only_resized_triangle(i_partner)%triangle, &
                this%short_interactions%visitable_cells_memento, this%short_interactions%&
                cells(ij_boxes(i_partner)))
            call short_interactions_reset(this%short_interactions%cells(ij_boxes(i_partner))%&
                neighbour_cells, only_resized_triangle(i_partner)%triangle, this%&
                short_interactions%cells(ij_boxes(i_partner))%visitable_cells)
            call this%dipolar_interactions_facades(ij_boxes(i_partner))%&
                reset(reset_real_pair(i_partner))
            call observables_energies_create(new_energies(i_partner), &
                size(this%mixture%components, 1))
        end do
    end subroutine Concrete_try

    subroutine Concrete_metropolis_algorithm(this, success, new_energies, ij_boxes, &
        boxes_size_ratio, energies)
        class(Boxes_Volume_Exchange), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Observables_Energies), intent(inout) :: new_energies(:)
        integer, intent(in) :: ij_boxes(:)
        real(DP), intent(in) :: boxes_size_ratio(:, :)
        type(Concrete_Observables_Energies), intent(in) :: energies(:)
    end subroutine Concrete_metropolis_algorithm

end module classes_boxes_volume_exchange
