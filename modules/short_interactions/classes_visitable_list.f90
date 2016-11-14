module classes_visitable_list

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_box_size, only: box_size_max_distance => max_distance
use types_particle_wrapper, only: Concrete_Particle
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_hard_contact, only: Abstract_Hard_Contact
use classes_pair_potential, only: Abstract_Pair_Potential
use classes_dipoles_neighbourhood, only: Abstract_Dipolar_Neighbourhood
use procedures_visit_condition, only: abstract_visit_condition
use module_list_node, only: Concrete_Linkable_Node, deallocate_list, increase_nodes_size

implicit none

private

    type, abstract, public :: Abstract_Visitable_List
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Component_Coordinates), pointer :: positions => null()
        class(Abstract_Component_Coordinates), pointer :: orientations => null()
        class(Abstract_Hard_Contact), pointer :: hard_contact => null()
        class(Abstract_Pair_Potential), pointer :: pair_potential => null()
        class(Abstract_Dipolar_Neighbourhood), pointer :: dipolar_neighbourhood => null()
        type(Concrete_Linkable_Node), pointer :: beginning => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: set => Abstract_set
        procedure :: add => Abstract_add
        procedure :: remove => Abstract_remove
        procedure :: visit_energy => Abstract_visit_energy
        procedure :: visit_contacts => Abstract_visit_contacts
        procedure :: visit_min_distance => Abstract_visit_min_distance
        procedure :: visit_dipolar_neighbours => Abstract_visit_dipolar_neighbours
    end type Abstract_Visitable_List

    type, extends(Abstract_Visitable_List), public :: Concrete_Visitable_List

    end type Concrete_Visitable_List

    type, extends(Abstract_Visitable_List), public :: Concrete_Visitable_Array
    private
        integer, allocatable :: nodes(:)
        integer :: num_nodes = 0
    contains
        procedure :: construct =>  Array_construct
        procedure :: destroy => Array_destroy
        procedure :: set => Array_set
        procedure :: add => Array_add
        procedure :: remove => Array_remove
        procedure :: visit_energy => Array_visit_energy
        procedure :: visit_contacts => Array_visit_contacts
        procedure :: visit_min_distance => Array_visit_min_distance
        procedure :: visit_dipolar_neighbours => Array_visit_dipolar_neighbours
    end type Concrete_Visitable_Array

    type, extends(Abstract_Visitable_List), public :: Null_Visitable_List
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: set => Null_set
        procedure :: add => Null_add
        procedure :: remove => Null_remove
        procedure :: visit_energy => Null_visit_energy
        procedure :: visit_contacts => Null_visit_contacts
        procedure :: visit_min_distance => Null_visit_min_distance
        procedure :: visit_dipolar_neighbours => Null_visit_dipolar_neighbours
    end type Null_Visitable_List

contains

!implementation Abstract_Visitable_List

    subroutine Abstract_construct(this, periodic_box, positions, orientations, hard_contact, &
        pair_potential, dipolar_neighbourhood)
        class(Abstract_Visitable_List), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: positions, orientations
        class(Abstract_Hard_Contact), target, intent(in) :: hard_contact
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential
        class(Abstract_Dipolar_Neighbourhood), target, intent(in) :: dipolar_neighbourhood

        type(Concrete_Linkable_Node), pointer :: current => null(), next => null()

        this%periodic_box => periodic_box
        this%positions => positions
        this%orientations => orientations
        this%hard_contact => hard_contact
        this%pair_potential => pair_potential
        this%dipolar_neighbourhood => dipolar_neighbourhood

        allocate(this%beginning)
        current => this%beginning
        current%i = 0
        allocate(next)
        next%i = 0
        current%next => next
        current => next
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Visitable_List), intent(inout) :: this

        call deallocate_list(this%beginning)
        this%dipolar_neighbourhood => null()
        this%pair_potential => null()
        this%hard_contact => null()
        this%orientations => null()
        this%positions => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_set(this, i_target, i_source)
        class(Abstract_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_target, i_source

        type(Concrete_Linkable_Node), pointer :: previous => null(), current => null(), &
            next => null()

        previous => this%beginning
        current => previous%next
        do
            next => current%next
            if (current%i == i_target) then
                current%i = i_source
                return
            else
                previous => current
            end if
            current => next
        end do
    end subroutine Abstract_set

    subroutine Abstract_add(this, i_particle)
        class(Abstract_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_particle

        type(Concrete_Linkable_Node), pointer :: previous => null(), new => null(), next => null()

        previous => this%beginning
        next => previous%next
        allocate(new)
        new%next => previous%next
        previous%next => new
        new%i = i_particle
    end subroutine Abstract_add

    subroutine Abstract_remove(this, i_particle)
        class(Abstract_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_particle

        type(Concrete_Linkable_Node), pointer :: previous => null(), current => null(), &
            next => null()

        previous => this%beginning
        current => previous%next
        do
            next => current%next
            if (current%i == i_particle) then
                previous%next => current%next
                deallocate(current)
                current => next
                return
            else
                previous => current
            end if
            current => next
        end do
    end subroutine Abstract_remove

    subroutine Abstract_visit_energy(this, overlap, energy, particle, visit_condition, i_exclude)
        class(Abstract_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        type(Concrete_Linkable_Node), pointer :: current => null(), next => null()
        real(DP) :: energy_i

        overlap = .false.
        energy = 0._DP
        current => this%beginning%next
        if (.not.associated(current%next)) return
        do
            next => current%next
            if (visit_condition(current%i, i_exclude)) then
                call this%pair_potential%meet(overlap, energy_i, this%periodic_box%&
                    distance(this%positions%get(current%i), particle%position))
                if (overlap) return
                energy = energy + energy_i
            end if
            if (.not.associated(next%next)) return
            current => next
        end do
    end subroutine Abstract_visit_energy

    subroutine Abstract_visit_contacts(this, overlap, contacts, particle, visit_condition, &
        i_exclude)
        class(Abstract_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: contacts
        type(Concrete_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        type(Concrete_Linkable_Node), pointer :: current => null(), next => null()
        real(DP) :: contact_i
        real(DP) :: min_distance

        overlap = .false.
        contacts = 0._DP
        current => this%beginning%next
        if (.not.associated(current%next)) return
        min_distance = this%pair_potential%get_min_distance()
        do
            next => current%next
            if (visit_condition(current%i, i_exclude)) then
                call this%hard_contact%meet(overlap, contact_i, min_distance, this%periodic_box%&
                    vector(this%positions%get(current%i), particle%position))
                if (overlap) return
                contacts = contacts + contact_i
            end if
            if (.not.associated(next%next)) return
            current => next
        end do
    end subroutine Abstract_visit_contacts

    subroutine Abstract_visit_min_distance(this, overlap, ratio, particle, visit_condition, &
        i_exclude)
        class(Abstract_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: ratio
        type(Concrete_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        type(Concrete_Linkable_Node), pointer :: current => null(), next => null()
        real(DP) :: min_distance
        logical :: can_overlap
        real(DP) :: ratio_i

        overlap = .false.
        min_distance = this%pair_potential%get_min_distance()
        ratio = box_size_max_distance(this%periodic_box%get_size()) / min_distance
        current => this%beginning%next
        if (.not.associated(current%next)) return
        do
            next => current%next
            if (visit_condition(current%i, i_exclude)) then
                call this%hard_contact%meet(can_overlap, overlap, ratio_i, min_distance, this%&
                    periodic_box%vector(this%positions%get(current%i), particle%position))
                if (can_overlap) then
                    if (overlap) return
                    if (ratio_i < ratio) ratio = ratio_i
                end if
            end if
            if (.not.associated(next%next)) return
            current => next
        end do
    end subroutine Abstract_visit_min_distance

    subroutine Abstract_visit_dipolar_neighbours(this, overlap, adjacency_matrix, particle, &
        visit_condition, i_exclude)
        class(Abstract_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        logical, intent(inout) :: adjacency_matrix(:, :)
        type(Concrete_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        type(Concrete_Linkable_Node), pointer :: current => null(), next => null()
        real(DP) :: min_distance

        overlap = .false.
        current => this%beginning%next
        if (.not.associated(current%next)) return
        min_distance = this%pair_potential%get_min_distance()
        do
            next => current%next
            if (visit_condition(current%i, i_exclude)) then
                call this%dipolar_neighbourhood%meet(overlap, &
                    adjacency_matrix(current%i, particle%i), min_distance, &
                    this%periodic_box%vector(this%positions%get(current%i), particle%position), &
                            this%orientations%get(current%i), particle%orientation)
                if (overlap) return
            end if
            if (.not.associated(next%next)) return
            current => next
        end do
    end subroutine Abstract_visit_dipolar_neighbours

!end implementation Abstract_Visitable_List

!implementation Concrete_Visitable_Array

    subroutine Array_construct(this, periodic_box, positions, orientations, hard_contact, &
        pair_potential, dipolar_neighbourhood)
        class(Concrete_Visitable_Array), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: positions, orientations
        class(Abstract_Hard_Contact), target, intent(in) :: hard_contact
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential
        class(Abstract_Dipolar_Neighbourhood), target, intent(in) :: dipolar_neighbourhood

        integer :: initial_size

        this%periodic_box => periodic_box
        this%positions => positions
        this%orientations => orientations
        this%hard_contact => hard_contact
        this%pair_potential => pair_potential
        this%dipolar_neighbourhood => dipolar_neighbourhood

        this%num_nodes = 0
        initial_size = 1
        allocate(this%nodes(initial_size))
    end subroutine Array_construct

    subroutine Array_destroy(this)
        class(Concrete_Visitable_Array), intent(inout) :: this

        if (allocated(this%nodes)) deallocate(this%nodes)
        this%dipolar_neighbourhood => null()
        this%pair_potential => null()
        this%hard_contact => null()
        this%orientations => null()
        this%positions => null()
        this%periodic_box => null()
    end subroutine Array_destroy

    subroutine Array_set(this, i_target, i_source)
        class(Concrete_Visitable_Array), intent(inout) :: this
        integer, intent(in) :: i_target, i_source

        integer :: i_node

        do i_node = 1, this%num_nodes
            if (this%nodes(i_node) == i_target) then
                this%nodes(i_node) = i_source
                return
            end if
        end do
    end subroutine Array_set

    subroutine Array_add(this, i_particle)
        class(Concrete_Visitable_Array), intent(inout) :: this
        integer, intent(in) :: i_particle

        this%num_nodes = this%num_nodes + 1
        if (size(this%nodes) < this%num_nodes) then
            call increase_nodes_size(this%nodes)
        end if
        this%nodes(this%num_nodes) = i_particle
    end subroutine Array_add

    subroutine Array_remove(this, i_particle)
        class(Concrete_Visitable_Array), intent(inout) :: this
        integer, intent(in) :: i_particle

        integer :: i_last

        if (i_particle /= this%nodes(this%num_nodes)) then
            i_last = this%nodes(this%num_nodes)
            call this%set(i_particle, i_last)
        end if
        this%num_nodes = this%num_nodes - 1
    end subroutine Array_remove

    subroutine Array_visit_energy(this, overlap, energy, particle, visit_condition, i_exclude)
        class(Concrete_Visitable_Array), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        real(DP) :: energy_i
        integer :: i_node

        overlap = .false.
        energy = 0._DP
        do i_node = 1, this%num_nodes
            if (.not.visit_condition(this%nodes(i_node), i_exclude)) cycle
            call this%pair_potential%meet(overlap, energy_i, this%periodic_box%&
                distance(this%positions%get(this%nodes(i_node)), particle%position))
            if (overlap) return
            energy = energy + energy_i
        end do
    end subroutine Array_visit_energy

    subroutine Array_visit_contacts(this, overlap, contacts, particle, visit_condition, i_exclude)
        class(Concrete_Visitable_Array), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: contacts
        type(Concrete_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        real(DP) :: contact_i
        real(DP) :: min_distance
        integer :: i_node

        min_distance = this%pair_potential%get_min_distance()
        overlap = .false.
        contacts = 0._DP
        do i_node = 1, this%num_nodes
            if (.not.visit_condition(this%nodes(i_node), i_exclude)) cycle
            call this%hard_contact%meet(overlap, contact_i, min_distance, this%periodic_box%&
                vector(this%positions%get(this%nodes(i_node)), particle%position))
            if (overlap) return
            contacts = contacts + contact_i
        end do
    end subroutine Array_visit_contacts

    subroutine Array_visit_min_distance(this, overlap, ratio, particle, visit_condition, i_exclude)
        class(Concrete_Visitable_Array), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: ratio
        type(Concrete_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        logical :: can_overlap
        real(DP) :: ratio_i, min_distance
        integer :: i_node

        overlap = .false.
        min_distance = this%pair_potential%get_min_distance()
        ratio = box_size_max_distance(this%periodic_box%get_size()) / min_distance
        do i_node = 1, this%num_nodes
            if (.not.visit_condition(this%nodes(i_node), i_exclude)) cycle
            call this%hard_contact%meet(can_overlap, overlap, ratio_i, min_distance, this%&
                periodic_box%vector(this%positions%get(this%nodes(i_node)), particle%position))
            if (.not. can_overlap) cycle
            if (overlap) return
            if (ratio_i < ratio) ratio = ratio_i
        end do
    end subroutine Array_visit_min_distance

    subroutine Array_visit_dipolar_neighbours(this, overlap, adjacency_matrix, particle, &
        visit_condition, i_exclude)
        class(Concrete_Visitable_Array), intent(in) :: this
        logical, intent(out) :: overlap
        logical, intent(inout) :: adjacency_matrix(:, :)
        type(Concrete_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        real(DP) :: min_distance
        integer :: i_node

        overlap = .false.
        min_distance = this%pair_potential%get_min_distance()
        do i_node = 1, this%num_nodes
            if (.not.visit_condition(this%nodes(i_node), i_exclude)) cycle
            call this%dipolar_neighbourhood%meet(overlap, &
                adjacency_matrix(this%nodes(i_node), particle%i), min_distance, &
                this%periodic_box%&
                    vector(this%positions%get(this%nodes(i_node)), particle%position), &
                    this%orientations%get(this%nodes(i_node)), particle%orientation)
            if (overlap) return
        end do
    end subroutine Array_visit_dipolar_neighbours

!end implementation Concrete_Visitable_Array

!implementation Null_Visitable_List

    subroutine Null_construct(this, periodic_box, positions, orientations, hard_contact, &
        pair_potential, dipolar_neighbourhood)
        class(Null_Visitable_List), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: positions, orientations
        class(Abstract_Hard_Contact), target, intent(in) :: hard_contact
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential
        class(Abstract_Dipolar_Neighbourhood), target, intent(in) :: dipolar_neighbourhood
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Visitable_List), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_set(this, i_target, i_source)
        class(Null_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_target, i_source
    end subroutine Null_set

    subroutine Null_visit_energy(this, overlap, energy, particle, visit_condition, i_exclude)
        class(Null_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude
        overlap = .false.; energy = 0._DP
    end subroutine Null_visit_energy

    subroutine Null_visit_contacts(this, overlap, contacts, particle, visit_condition, i_exclude)
        class(Null_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: contacts
        type(Concrete_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude
        overlap = .false.; contacts = 0._DP
    end subroutine Null_visit_contacts

    subroutine Null_visit_min_distance(this, overlap, ratio, particle, visit_condition, i_exclude)
        class(Null_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: ratio
        type(Concrete_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude
        overlap = .false.; ratio = 0._DP
    end subroutine Null_visit_min_distance

    subroutine Null_visit_dipolar_neighbours(this, overlap, adjacency_matrix, particle, &
        visit_condition, i_exclude)
        class(Null_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        logical, intent(inout) :: adjacency_matrix(:, :)
        type(Concrete_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude
        overlap = .false.
    end subroutine Null_visit_dipolar_neighbours

    subroutine Null_add(this, i_particle)
        class(Null_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_add

    subroutine Null_remove(this, i_particle)
        class(Null_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_remove

!end implementation Null_Visitable_List

end module classes_visitable_list
