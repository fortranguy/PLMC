module classes_visitable_cells

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use data_cells, only: nums_local_cells
use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_component_coordinates, only: Abstract_Component_Coordinates
use types_temporary_particle, only: Concrete_Temporary_Particle
use classes_hard_contact, only: Abstract_Hard_Contact
use procedures_visit_condition, only: abstract_visit_condition
use classes_visitable_list, only: Abstract_Visitable_List
use classes_pair_potential, only: Abstract_Pair_Potential
use classes_neighbour_cells, only: Abstract_Neighbour_Cells

implicit none

private

    type, abstract, public :: Abstract_Visitable_Cells
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Component_Coordinates), pointer :: positions => null()
        class(Abstract_Hard_Contact), pointer :: hard_contact => null()
        class(Abstract_Pair_Potential), pointer :: pair_potential => null()
        class(Abstract_Neighbour_Cells), pointer :: neighbour_cells => null()
        class(Abstract_Visitable_List), allocatable :: visitable_lists(:, :, :), list_mold
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset => Abstract_reset
        procedure :: set => Abstract_set
        procedure :: visit_energy => Abstract_visit_energy
        procedure :: visit_contacts => Abstract_visit_contacts
        procedure :: visit_min_distance => Abstract_visit_min_distance
        procedure :: translate => Abstract_translate
        procedure :: add => Abstract_add
        procedure :: remove => Abstract_remove
        procedure, private :: construct_visitable_lists => Abstract_construct_visitable_lists
        procedure, private :: destroy_visitable_lists => Abstract_destroy_visitable_lists
        procedure, private :: fill_with_particles => Abstract_fill_with_particles
    end type Abstract_Visitable_Cells

    type, extends(Abstract_Visitable_Cells), public :: Concrete_Visitable_Cells

    end type Concrete_Visitable_Cells

    type, extends(Abstract_Visitable_Cells), public :: Null_Visitable_Cells
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: reset => Null_reset
        procedure :: set => Null_set
        procedure :: visit_energy => Null_visit_energy
        procedure :: visit_contacts => Null_visit_contacts
        procedure :: visit_min_distance => Null_visit_min_distance
        procedure :: translate => Null_translate
        procedure :: add => Null_add
        procedure :: remove => Null_remove
    end type Null_Visitable_Cells

contains

!implementation Abstract_Visitable_Cells

    subroutine Abstract_construct(this, periodic_box, positions, hard_contact, pair_potential, &
        list_mold)
        class(Abstract_Visitable_Cells), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: positions
        class(Abstract_Hard_Contact), target, intent(in) :: hard_contact
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential
        class(Abstract_Visitable_List), intent(in) :: list_mold

        this%periodic_box => periodic_box
        this%positions => positions
        this%hard_contact => hard_contact
        this%pair_potential => pair_potential
        allocate(this%list_mold, mold=list_mold)
    end subroutine Abstract_construct

    subroutine Abstract_reset(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        call this%destroy_visitable_lists()
        call this%construct_visitable_lists()
        call this%fill_with_particles()
    end subroutine Abstract_reset

    subroutine Abstract_set(this, neighbour_cells)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        class(Abstract_Neighbour_Cells), target, intent(in) :: neighbour_cells

        this%neighbour_cells => neighbour_cells
    end subroutine Abstract_set

    subroutine Abstract_construct_visitable_lists(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        integer, dimension(num_dimensions) :: global_lbounds, global_ubounds
        integer :: global_i1, global_i2, global_i3

        global_lbounds = this%neighbour_cells%get_global_lbounds()
        global_ubounds = this%neighbour_cells%get_global_ubounds()
        allocate(this%visitable_lists(global_lbounds(1):global_ubounds(1), &
                                      global_lbounds(2):global_ubounds(2), &
                                      global_lbounds(3):global_ubounds(3)), &
                                      mold=this%list_mold)
        do global_i3 = lbound(this%visitable_lists, 3), ubound(this%visitable_lists, 3)
        do global_i2 = lbound(this%visitable_lists, 2), ubound(this%visitable_lists, 2)
        do global_i1 = lbound(this%visitable_lists, 1), ubound(this%visitable_lists, 1)
            call this%visitable_lists(global_i1, global_i2, global_i3)%construct(this%periodic_box,&
                this%positions, this%hard_contact)
        end do
        end do
        end do
    end subroutine Abstract_construct_visitable_lists

    subroutine Abstract_fill_with_particles(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        type(Concrete_Temporary_Particle) :: particle
        integer :: i_particle

        do i_particle = 1, this%positions%get_num()
            particle%i = i_particle
            particle%position = this%positions%get(particle%i)
            call this%add(particle)
        end do
    end subroutine Abstract_fill_with_particles

    subroutine Abstract_destroy(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        call this%destroy_visitable_lists()
        if (allocated(this%list_mold)) deallocate(this%list_mold)
        this%neighbour_cells => null()
        this%pair_potential => null()
        this%hard_contact => null()
        this%positions => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_destroy_visitable_lists(this)
        class(Abstract_Visitable_Cells), intent(inout) :: this

        integer :: global_i1, global_i2, global_i3

        if (.not. allocated(this%visitable_lists)) return
        do global_i3 = lbound(this%visitable_lists, 3), ubound(this%visitable_lists, 3)
        do global_i2 = lbound(this%visitable_lists, 2), ubound(this%visitable_lists, 2)
        do global_i1 = lbound(this%visitable_lists, 1), ubound(this%visitable_lists, 1)
            call this%visitable_lists(global_i1, global_i2, global_i3)%destroy()
        end do
        end do
        end do
        deallocate(this%visitable_lists)
    end subroutine Abstract_destroy_visitable_lists

    subroutine Abstract_visit_energy(this, overlap, energy, particle, visit_condition, i_exclude)
        class(Abstract_Visitable_Cells), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        real(DP) :: energy_i
        integer, dimension(num_dimensions) :: ijk_cell, ijk_local_cell
        logical :: at_bottom_layer, at_top_layer
        integer :: local_i1, local_i2, local_i3

        ijk_cell = this%neighbour_cells%index(particle%position)
        at_bottom_layer = (ijk_cell(3) == lbound(this%visitable_lists, 3))
        at_top_layer = (ijk_cell(3) == ubound(this%visitable_lists, 3))
        energy = 0._DP
        do local_i3 = -nums_local_cells(3)/2, nums_local_cells(3)/2
            if (this%neighbour_cells%skip(at_bottom_layer, at_top_layer, local_i3)) cycle
        do local_i2 = -nums_local_cells(2)/2, nums_local_cells(2)/2
        do local_i1 = -nums_local_cells(1)/2, nums_local_cells(1)/2
            ijk_local_cell = this%neighbour_cells%get(local_i1, local_i2, local_i3, ijk_cell(1), &
                ijk_cell(2), ijk_cell(3))
            call this%visitable_lists(ijk_local_cell(1), ijk_local_cell(2), ijk_local_cell(3))%&
                visit(overlap, energy_i, particle, this%pair_potential, visit_condition, i_exclude)
            if (overlap) return
            energy = energy + energy_i
        end do
        end do
        end do
    end subroutine Abstract_visit_energy

    subroutine Abstract_visit_contacts(this, overlap, contacts, particle, visit_condition, &
        i_exclude)
        class(Abstract_Visitable_Cells), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: contacts
        type(Concrete_Temporary_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        real(DP) :: contacts_i, min_distance
        integer, dimension(num_dimensions) :: ijk_cell, ijk_local_cell
        logical :: at_bottom_layer, at_top_layer
        integer :: local_i1, local_i2, local_i3

        min_distance = this%pair_potential%get_min_distance()
        ijk_cell = this%neighbour_cells%index(particle%position)
        at_bottom_layer = (ijk_cell(3) == lbound(this%visitable_lists, 3))
        at_top_layer = (ijk_cell(3) == ubound(this%visitable_lists, 3))
        contacts = 0._DP
        do local_i3 = -nums_local_cells(3)/2, nums_local_cells(3)/2
            if (this%neighbour_cells%skip(at_bottom_layer, at_top_layer, local_i3)) cycle
        do local_i2 = -nums_local_cells(2)/2, nums_local_cells(2)/2
        do local_i1 = -nums_local_cells(1)/2, nums_local_cells(1)/2
            ijk_local_cell = this%neighbour_cells%get(local_i1, local_i2, local_i3, ijk_cell(1), &
                ijk_cell(2), ijk_cell(3))
            call this%visitable_lists(ijk_local_cell(1), ijk_local_cell(2), ijk_local_cell(3))%&
                visit(overlap, contacts_i, particle, min_distance, visit_condition, i_exclude)
            if (overlap) return
            contacts = contacts + contacts_i
        end do
        end do
        end do
    end subroutine Abstract_visit_contacts

    subroutine Abstract_visit_min_distance(this, can_overlap, overlap, ratio, particle, &
        visit_condition, i_exclude)
        class(Abstract_Visitable_Cells), intent(in) :: this
        logical, intent(out) :: can_overlap, overlap
        real(DP), intent(out) :: ratio
        type(Concrete_Temporary_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        real(DP) :: ratio_i, min_distance, max_distance
        integer, dimension(num_dimensions) :: ijk_cell, ijk_local_cell
        logical :: at_bottom_layer, at_top_layer
        integer :: local_i1, local_i2, local_i3

        min_distance = this%pair_potential%get_min_distance()
        max_distance = this%periodic_box%get_max_distance()
        ratio = max_distance / min_distance
        ijk_cell = this%neighbour_cells%index(particle%position)
        at_bottom_layer = (ijk_cell(3) == lbound(this%visitable_lists, 3))
        at_top_layer = (ijk_cell(3) == ubound(this%visitable_lists, 3))
        do local_i3 = -nums_local_cells(3)/2, nums_local_cells(3)/2
            if (this%neighbour_cells%skip(at_bottom_layer, at_top_layer, local_i3)) cycle
        do local_i2 = -nums_local_cells(2)/2, nums_local_cells(2)/2
        do local_i1 = -nums_local_cells(1)/2, nums_local_cells(1)/2
            ijk_local_cell = this%neighbour_cells%get(local_i1, local_i2, local_i3, ijk_cell(1), &
                ijk_cell(2), ijk_cell(3))
            call this%visitable_lists(ijk_local_cell(1), ijk_local_cell(2), ijk_local_cell(3))%&
                visit(can_overlap, overlap, ratio_i, particle, min_distance, max_distance, &
                    visit_condition, i_exclude)
            if (.not. can_overlap) cycle
            if (overlap) return
            if (ratio_i < ratio) ratio = ratio_i
        end do
        end do
        end do
    end subroutine Abstract_visit_min_distance

    !> No check: to_position & from%position are assumed to be within accessible_domain
    subroutine Abstract_translate(this, to_position, from)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        real(DP), intent(in) :: to_position(:)
        type(Concrete_Temporary_Particle), intent(in) :: from

        integer, dimension(num_dimensions) :: from_ijk_cell, to_ijk_cell

        from_ijk_cell = this%neighbour_cells%index(from%position)
        to_ijk_cell = this%neighbour_cells%index(to_position)
        if (any(from_ijk_cell /= to_ijk_cell)) then
            call this%visitable_lists(from_ijk_cell(1), from_ijk_cell(2), from_ijk_cell(3))%&
                remove(from%i)
            call this%visitable_lists(to_ijk_cell(1), to_ijk_cell(2), to_ijk_cell(3))%add(from%i)
        end if
    end subroutine Abstract_translate

    subroutine Abstract_add(this, particle)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle

        integer :: ijk_cell(num_dimensions)

        if (.not.this%neighbour_cells%is_inside(particle%position)) then
            call error_exit("Abstract_Visitable_Cells: add: particle%position is outside "//&
                "accessible_domain.")
        end if
        ijk_cell = this%neighbour_cells%index(particle%position)
        call this%visitable_lists(ijk_cell(1), ijk_cell(2), ijk_cell(3))%add(particle%i)
    end subroutine Abstract_add

    subroutine Abstract_remove(this, particle)
        class(Abstract_Visitable_Cells), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle

        integer :: ijk_cell(num_dimensions)

        ijk_cell = this%neighbour_cells%index(particle%position)
        call this%visitable_lists(ijk_cell(1), ijk_cell(2), ijk_cell(3))%remove(particle%i)
        if (particle%i < this%positions%get_num()) then
            ijk_cell = this%neighbour_cells%index(this%positions%get(this%positions%&
                get_num()))
            call this%visitable_lists(ijk_cell(1), ijk_cell(2), ijk_cell(3))%set(this%positions%&
                get_num(), particle%i)
        end if
    end subroutine Abstract_remove

!end implementation Abstract_Visitable_Cells

!implementation Null_Visitable_Cells

    subroutine Null_construct(this, periodic_box, positions, hard_contact, pair_potential, &
        list_mold)
        class(Null_Visitable_Cells), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: positions
        class(Abstract_Hard_Contact), target, intent(in) :: hard_contact
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential
        class(Abstract_Visitable_List), intent(in) :: list_mold
    end subroutine Null_construct

    subroutine Null_reset(this)
        class(Null_Visitable_Cells), intent(inout) :: this
    end subroutine Null_reset

    subroutine Null_set(this, neighbour_cells)
        class(Null_Visitable_Cells), intent(inout) :: this
        class(Abstract_Neighbour_Cells), target, intent(in) :: neighbour_cells
    end subroutine Null_set

    subroutine Null_destroy(this)
        class(Null_Visitable_Cells), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_visit_energy(this, overlap, energy, particle, visit_condition, i_exclude)
        class(Null_Visitable_Cells), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude
        overlap = .false.
        energy = 0._DP
    end subroutine Null_visit_energy

    subroutine Null_visit_contacts(this, overlap, contacts, particle, visit_condition, i_exclude)
        class(Null_Visitable_Cells), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: contacts
        type(Concrete_Temporary_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude
        overlap = .false.
        contacts = 0._DP
    end subroutine Null_visit_contacts

    subroutine Null_visit_min_distance(this, can_overlap, overlap, ratio, particle, &
        visit_condition, i_exclude)
        class(Null_Visitable_Cells), intent(in) :: this
        logical, intent(out) :: can_overlap, overlap
        real(DP), intent(out) :: ratio
        type(Concrete_Temporary_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude
        can_overlap = .false.; overlap = .false.
        ratio = 0._DP
    end subroutine Null_visit_min_distance

    subroutine Null_translate(this, to_position, from)
        class(Null_Visitable_Cells), intent(inout) :: this
        real(DP), intent(in) :: to_position(:)
        type(Concrete_Temporary_Particle), intent(in) :: from
    end subroutine Null_translate

    subroutine Null_add(this, particle)
        class(Null_Visitable_Cells), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle
    end subroutine Null_add

    subroutine Null_remove(this, particle)
        class(Null_Visitable_Cells), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: particle
    end subroutine Null_remove

!end implementation Null_Visitable_Cells

end module classes_visitable_cells
