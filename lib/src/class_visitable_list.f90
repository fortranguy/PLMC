module class_visitable_list

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_temporary_particle, only: Concrete_Temporary_Particle
use module_list_node, only: Concrete_Linkable_Node, &
    deallocate_list, increase_nodes_size
use class_periodic_box, only: Abstract_Periodic_Box
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Visitable_List
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box
        class(Abstract_Component_Coordinates), pointer :: component_positions
        type(Concrete_Linkable_Node), pointer :: beginning
    contains
        procedure :: construct => Abstract_Visitable_List_construct
        procedure :: destroy => Abstract_Visitable_List_destroy
        procedure :: set => Abstract_Visitable_List_set
        procedure :: visit => Abstract_Visitable_List_visit
        procedure :: add => Abstract_Visitable_List_add
        procedure :: remove => Abstract_Visitable_List_remove
    end type Abstract_Visitable_List

    type, extends(Abstract_Visitable_List), public :: Concrete_Visitable_List

    end type Concrete_Visitable_List

    type, extends(Abstract_Visitable_List), public :: Concrete_Visitable_Array
    private
        integer, allocatable :: nodes(:)
        integer :: num_nodes
    contains
        procedure :: construct =>  Concrete_Visitable_Array_construct
        procedure :: destroy => Concrete_Visitable_Array_destroy
        procedure :: set => Concrete_Visitable_Array_set
        procedure :: visit => Concrete_Visitable_Array_visit
        procedure :: add => Concrete_Visitable_Array_add
        procedure :: remove => Concrete_Visitable_Array_remove
    end type Concrete_Visitable_Array

    type, extends(Abstract_Visitable_List), public :: Null_Visitable_List
    contains
        procedure :: construct => Null_Visitable_List_construct
        procedure :: destroy => Null_Visitable_List_destroy
        procedure :: set => Null_Visitable_List_set
        procedure :: visit => Null_Visitable_List_visit
        procedure :: add => Null_Visitable_List_add
        procedure :: remove => Null_Visitable_List_remove
    end type Null_Visitable_List

contains

!implementation Abstract_Visitable_List

    subroutine Abstract_Visitable_List_construct(this, periodic_box, component_positions)
        class(Abstract_Visitable_List), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: component_positions

        type(Concrete_Linkable_Node), pointer :: current => null(), next => null()

        this%periodic_box => periodic_box
        this%component_positions => component_positions

        allocate(this%beginning)
        current => this%beginning
        current%i = 0
        allocate(next)
        next%i = 0
        current%next => next
        current => next
    end subroutine Abstract_Visitable_List_construct

    subroutine Abstract_Visitable_List_destroy(this)
        class(Abstract_Visitable_List), intent(inout) :: this

        call deallocate_list(this%beginning)
        this%component_positions => null()
        this%periodic_box => null()
    end subroutine Abstract_Visitable_List_destroy

    subroutine Abstract_Visitable_List_set(this, i_target, i_particle)
        class(Abstract_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_target, i_particle

        type(Concrete_Linkable_Node), pointer :: previous => null(), current => null(), &
            next => null()

        previous => this%beginning
        current => previous%next
        do
            next => current%next
            if (current%i == i_target) then
                current%i = i_particle
                return
            else
                previous => current
            end if
            current => next
        end do
    end subroutine Abstract_Visitable_List_set

    subroutine Abstract_Visitable_List_visit(this, overlap, energy, particle, pair_potential, &
        i_exclude)
        class(Abstract_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        integer, intent(in) :: i_exclude

        type(Concrete_Linkable_Node), pointer :: current => null(), next => null()
        real(DP) :: energy_i, distance

        overlap = .false.
        energy = 0._DP
        current => this%beginning%next
        if (.not.associated(current%next)) return
        do
            next => current%next
            if (current%i /= i_exclude) then
                distance = this%periodic_box%distance(particle%position, &
                    this%component_positions%get(current%i))
                call pair_potential%meet(overlap, energy_i, distance)
                if (overlap) return
                energy = energy + energy_i
            end if
            if (.not.associated(next%next)) return
            current => next
        end do
    end subroutine Abstract_Visitable_List_visit

    subroutine Abstract_Visitable_List_add(this, i_particle)
        class(Abstract_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_particle

        type(Concrete_Linkable_Node), pointer :: previous => null(), new => null(), next => null()

        previous => this%beginning
        next => previous%next
        allocate(new)
        new%next => previous%next
        previous%next => new
        new%i = i_particle
    end subroutine Abstract_Visitable_List_add

    subroutine Abstract_Visitable_List_remove(this, i_particle)
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
    end subroutine Abstract_Visitable_List_remove

!end implementation Abstract_Visitable_List

!implementation Concrete_Visitable_Array

    subroutine Concrete_Visitable_Array_construct(this, periodic_box, component_positions)
        class(Concrete_Visitable_Array), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: component_positions

        integer :: initial_size

        this%periodic_box => periodic_box
        this%component_positions => component_positions
        this%num_nodes = 0
        initial_size = 1
        allocate(this%nodes(initial_size))
    end subroutine Concrete_Visitable_Array_construct

    subroutine Concrete_Visitable_Array_destroy(this)
        class(Concrete_Visitable_Array), intent(inout) :: this

        if (allocated(this%nodes)) deallocate(this%nodes)
        this%component_positions => null()
        this%periodic_box => null()
    end subroutine Concrete_Visitable_Array_destroy

    subroutine Concrete_Visitable_Array_set(this, i_target, i_particle)
        class(Concrete_Visitable_Array), intent(inout) :: this
        integer, intent(in) :: i_target, i_particle

        integer :: i_node

        do i_node = 1, this%num_nodes
            if (this%nodes(i_node) == i_target) then
                this%nodes(i_node) = i_particle
                return
            end if
        end do
    end subroutine Concrete_Visitable_Array_set

    subroutine Concrete_Visitable_Array_visit(this, overlap, energy, particle, pair_potential, &
        i_exclude)
        class(Concrete_Visitable_Array), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        integer, intent(in) :: i_exclude

        real(DP) :: energy_i, distance
        integer :: i_node

        overlap = .false.
        energy = 0._DP
        do i_node = 1, this%num_nodes
            if (this%nodes(i_node) == i_exclude) cycle
            distance = this%periodic_box%distance(particle%position, this%component_positions%&
                get(this%nodes(i_node)))
            call pair_potential%meet(overlap, energy_i, distance)
            if (overlap) return
            energy = energy + energy_i
        end do
    end subroutine Concrete_Visitable_Array_visit

    subroutine Concrete_Visitable_Array_add(this, i_particle)
        class(Concrete_Visitable_Array), intent(inout) :: this
        integer, intent(in) :: i_particle

        this%num_nodes = this%num_nodes + 1
        if (size(this%nodes) < this%num_nodes) then
            call increase_nodes_size(this%nodes)
        end if
        this%nodes(this%num_nodes) = i_particle
    end subroutine Concrete_Visitable_Array_add

    subroutine Concrete_Visitable_Array_remove(this, i_particle)
        class(Concrete_Visitable_Array), intent(inout) :: this
        integer, intent(in) :: i_particle

        integer :: i_last

        if (i_particle /= this%nodes(this%num_nodes)) then
            i_last = this%nodes(this%num_nodes)
            call this%set(i_particle, i_last)
        end if
        this%num_nodes = this%num_nodes - 1
    end subroutine Concrete_Visitable_Array_remove

!end implementation Concrete_Visitable_Array

!implementation Null_Visitable_List

    subroutine Null_Visitable_List_construct(this, periodic_box, component_positions)
        class(Null_Visitable_List), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: component_positions
    end subroutine Null_Visitable_List_construct

    subroutine Null_Visitable_List_destroy(this)
        class(Null_Visitable_List), intent(inout) :: this
    end subroutine Null_Visitable_List_destroy

    subroutine Null_Visitable_List_set(this, i_target, i_particle)
        class(Null_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_target, i_particle
    end subroutine Null_Visitable_List_set

    subroutine Null_Visitable_List_visit(this, overlap, energy, particle, pair_potential, i_exclude)
        class(Null_Visitable_List), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        integer, intent(in) :: i_exclude
        overlap = .false.
        energy = 0._DP
    end subroutine Null_Visitable_List_visit

    subroutine Null_Visitable_List_add(this, i_particle)
        class(Null_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_Visitable_List_add

    subroutine Null_Visitable_List_remove(this, i_particle)
        class(Null_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_Visitable_List_remove

!end implementation Null_Visitable_List

end module class_visitable_list
