module classes_visitable_list

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use classes_periodic_box, only: Abstract_Periodic_Box
use types_temporary_particle, only: Concrete_Temporary_Particle
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_hard_contact, only: Abstract_Hard_Contact
use classes_pair_potential, only: Abstract_Pair_Potential
use procedures_visit_condition, only: abstract_visit_condition
use module_list_node, only: Concrete_Linkable_Node, deallocate_list, increase_nodes_size

implicit none

private

    type, abstract, public :: Abstract_Visitable_List
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Component_Coordinates), pointer :: positions => null()
        class(Abstract_Hard_Contact), pointer :: hard_contact => null()
        type(Concrete_Linkable_Node), pointer :: beginning => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        generic :: visit => visit_energy, visit_contacts, visit_min_distance
        procedure :: set => Abstract_set
        procedure :: add => Abstract_add
        procedure :: remove => Abstract_remove
        procedure, private :: visit_energy => Abstract_visit_energy
        procedure, private :: visit_contacts => Abstract_visit_contacts
        procedure, private :: visit_min_distance => Abstract_visit_min_distance
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
        procedure, private :: visit_energy => Array_visit_energy
        procedure, private :: visit_contacts => Array_visit_contacts
        procedure, private :: visit_min_distance => Array_visit_min_distance
    end type Concrete_Visitable_Array

    type, extends(Abstract_Visitable_List), public :: Null_Visitable_List
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: set => Null_set
        procedure :: add => Null_add
        procedure :: remove => Null_remove
        procedure, private :: visit_energy => Null_visit_energy
        procedure, private :: visit_contacts => Null_visit_contacts
        procedure, private :: visit_min_distance => Null_visit_min_distance
    end type Null_Visitable_List

contains

!implementation Abstract_Visitable_List

    subroutine Abstract_construct(this, periodic_box, positions, hard_contact)
        class(Abstract_Visitable_List), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: positions
        class(Abstract_Hard_Contact), target, intent(in) :: hard_contact

        type(Concrete_Linkable_Node), pointer :: current => null(), next => null()

        this%periodic_box => periodic_box
        this%positions => positions
        this%hard_contact => hard_contact

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
        this%hard_contact => null()
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

    subroutine Abstract_visit_energy(this, overlap, energy, particle, pair_potential, &
        visit_condition, i_exclude)
        class(Abstract_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        type(Concrete_Linkable_Node), pointer :: current => null(), next => null()
        real(DP) :: energy_i, distance

        overlap = .false.
        energy = 0._DP
        current => this%beginning%next
        if (.not.associated(current%next)) return
        do
            next => current%next
            if (visit_condition(current%i, i_exclude)) then
                distance = this%periodic_box%distance(particle%position, this%positions%&
                    get(current%i))
                call pair_potential%meet(overlap, energy_i, distance)
                if (overlap) return
                energy = energy + energy_i
            end if
            if (.not.associated(next%next)) return
            current => next
        end do
    end subroutine Abstract_visit_energy

    subroutine Abstract_visit_contacts(this, overlap, contacts, particle, min_distance, &
        visit_condition, i_exclude)
        class(Abstract_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: contacts
        type(Concrete_Temporary_Particle), intent(in) :: particle
        real(DP), intent(in) :: min_distance
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        type(Concrete_Linkable_Node), pointer :: current => null(), next => null()
        real(DP) :: contact_i, vector(num_dimensions)

        overlap = .false.
        contacts = 0._DP
        current => this%beginning%next
        if (.not.associated(current%next)) return
        do
            next => current%next
            if (visit_condition(current%i, i_exclude)) then
                vector = this%periodic_box%vector(particle%position, this%positions%get(current%i))
                call this%hard_contact%meet(overlap, contact_i, min_distance, vector)
                if (overlap) return
                contacts = contacts + contact_i
            end if
            if (.not.associated(next%next)) return
            current => next
        end do
    end subroutine Abstract_visit_contacts

    subroutine Abstract_visit_min_distance(this, overlap, ratio, particle, min_distance, &
        max_distance, visit_condition, i_exclude)
        class(Abstract_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: ratio
        type(Concrete_Temporary_Particle), intent(in) :: particle
        real(DP), intent(in) :: min_distance, max_distance
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        real(DP) :: ratio_i
        real(DP) :: vector(num_dimensions)
        logical :: can_overlap
        type(Concrete_Linkable_Node), pointer :: current => null(), next => null()

        overlap = .false.
        ratio = max_distance / min_distance
        current => this%beginning%next
        if (.not.associated(current%next)) return
        do
            next => current%next
            if (visit_condition(current%i, i_exclude)) then
                vector = this%periodic_box%vector(particle%position, this%positions%get(current%i))
                call this%hard_contact%meet(can_overlap, overlap, ratio_i, min_distance, vector)
                if (can_overlap) then
                    if (overlap) return
                    if (ratio_i < ratio) ratio = ratio_i
                end if
            end if
            if (.not.associated(next%next)) return
            current => next
        end do
    end subroutine Abstract_visit_min_distance

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

!end implementation Abstract_Visitable_List

!implementation Concrete_Visitable_Array

    subroutine Array_construct(this, periodic_box, positions, hard_contact)
        class(Concrete_Visitable_Array), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: positions
        class(Abstract_Hard_Contact), target, intent(in) :: hard_contact

        integer :: initial_size

        this%periodic_box => periodic_box
        this%positions => positions
        this%hard_contact => hard_contact
        this%num_nodes = 0
        initial_size = 1
        allocate(this%nodes(initial_size))
    end subroutine Array_construct

    subroutine Array_destroy(this)
        class(Concrete_Visitable_Array), intent(inout) :: this

        if (allocated(this%nodes)) deallocate(this%nodes)
        this%hard_contact => null()
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

    subroutine Array_visit_energy(this, overlap, energy, particle, pair_potential, visit_condition,&
        i_exclude)
        class(Concrete_Visitable_Array), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        real(DP) :: energy_i, distance
        integer :: i_node

        overlap = .false.
        energy = 0._DP
        do i_node = 1, this%num_nodes
            if (.not.visit_condition(this%nodes(i_node), i_exclude)) cycle
            distance = this%periodic_box%distance(particle%position, this%positions%&
                get(this%nodes(i_node)))
            call pair_potential%meet(overlap, energy_i, distance)
            if (overlap) return
            energy = energy + energy_i
        end do
    end subroutine Array_visit_energy

    subroutine Array_visit_contacts(this, overlap, contacts, particle, min_distance, &
        visit_condition, i_exclude)
        class(Concrete_Visitable_Array), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: contacts
        type(Concrete_Temporary_Particle), intent(in) :: particle
        real(DP), intent(in) :: min_distance
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        real(DP) :: contact_i, vector(num_dimensions)
        integer :: i_node

        overlap = .false.
        contacts = 0._DP
        do i_node = 1, this%num_nodes
            if (.not.visit_condition(this%nodes(i_node), i_exclude)) cycle
            vector = this%periodic_box%vector(particle%position, this%positions%get(this%&
                nodes(i_node)))
            call this%hard_contact%meet(overlap, contact_i, min_distance, vector)
            if (overlap) return
            contacts = contacts + contact_i
        end do
    end subroutine Array_visit_contacts

    subroutine Array_visit_min_distance(this, overlap, ratio, particle, min_distance, max_distance,&
        visit_condition, i_exclude)
        class(Concrete_Visitable_Array), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: ratio
        type(Concrete_Temporary_Particle), intent(in) :: particle
        real(DP), intent(in) :: min_distance, max_distance
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        real(DP) :: ratio_i
        real(DP) :: vector(num_dimensions)
        logical :: can_overlap
        integer :: i_node

        overlap = .false.
        ratio = max_distance / min_distance
        do i_node = 1, this%num_nodes
            if (.not.visit_condition(this%nodes(i_node), i_exclude)) cycle
            vector = this%periodic_box%vector(particle%position, this%positions%get(this%&
                nodes(i_node)))
            call this%hard_contact%meet(can_overlap, overlap, ratio_i, min_distance, vector)
            if (.not. can_overlap) cycle
            if (overlap) return
            if (ratio_i < ratio) ratio = ratio_i
        end do
    end subroutine Array_visit_min_distance

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

!end implementation Concrete_Visitable_Array

!implementation Null_Visitable_List

    subroutine Null_construct(this, periodic_box, positions, hard_contact)
        class(Null_Visitable_List), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: positions
        class(Abstract_Hard_Contact), target, intent(in) :: hard_contact
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Visitable_List), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_set(this, i_target, i_source)
        class(Null_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_target, i_source
    end subroutine Null_set

    subroutine Null_visit_energy(this, overlap, energy, particle, pair_potential, visit_condition, &
        i_exclude)
        class(Null_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude
        overlap = .false.
        energy = 0._DP
    end subroutine Null_visit_energy

    subroutine Null_visit_contacts(this, overlap, contacts, particle, min_distance, &
        visit_condition, i_exclude)
        class(Null_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: contacts
        type(Concrete_Temporary_Particle), intent(in) :: particle
        real(DP), intent(in) :: min_distance
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude
        overlap = .false.
        contacts = 0._DP
    end subroutine Null_visit_contacts

    subroutine Null_visit_min_distance(this, overlap, ratio, particle, min_distance, &
        max_distance, visit_condition, i_exclude)
        class(Null_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: ratio
        type(Concrete_Temporary_Particle), intent(in) :: particle
        real(DP), intent(in) :: min_distance, max_distance
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude
        overlap = .false.
        ratio = 0._DP
    end subroutine Null_visit_min_distance

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
