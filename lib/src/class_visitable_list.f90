module class_visitable_list

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_particle, only: Concrete_Particle
use module_nodes, only: Concrete_Node, Concrete_Linkable_Node, &
    deallocate_list, increase_nodes_size
use class_periodic_box, only: Abstract_Periodic_Box
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Visitable_List
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box
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
        type(Concrete_Node), allocatable :: nodes(:)
        integer :: num_nodes
    contains
        procedure :: construct =>  Concrete_Visitable_Array_construct
        procedure :: destroy => Concrete_Visitable_Array_destroy
        procedure :: set => Concrete_Visitable_Array_set
        procedure :: visit => Concrete_Visitable_Array_visit
        procedure :: add => Concrete_Visitable_Array_add
        procedure :: remove => Concrete_Visitable_Array_remove
    end type Concrete_Visitable_Array

contains

!implementation Abstract_Visitable_List

    subroutine Abstract_Visitable_List_construct(this, periodic_box)
        class(Abstract_Visitable_List), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box

        type(Concrete_Linkable_Node), pointer :: current => null(), next => null()

        this%periodic_box => periodic_box

        allocate(this%beginning)
        current => this%beginning
        current%i = 0
        current%position = 0._DP
        allocate(next)
        next%i = 0
        next%position = 0._DP
        current%next => next
        current => next
    end subroutine Abstract_Visitable_List_construct

    subroutine Abstract_Visitable_List_destroy(this)
        class(Abstract_Visitable_List), intent(inout) :: this

        call deallocate_list(this%beginning)
        this%periodic_box => null()
    end subroutine Abstract_Visitable_List_destroy

    subroutine Abstract_Visitable_List_set(this, i_target, particle)
        class(Abstract_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_target
        type(Concrete_Particle), intent(in) :: particle

        type(Concrete_Linkable_Node), pointer :: previous => null(), current => null(), next => null()

        previous => this%beginning
        current => previous%next
        do
            next => current%next
            if (current%i == i_target) then
                current%i = particle%i
                current%position = particle%position
                return
            else
                previous => current
            end if
            current => next
        end do
    end subroutine Abstract_Visitable_List_set

    subroutine Abstract_Visitable_List_visit(this, particle, pair_potential, overlap, energy)
        class(Abstract_Visitable_List), intent(in) :: this
        type(Concrete_Particle), intent(in) :: particle
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy

        type(Concrete_Linkable_Node), pointer :: current => null(), next => null()
        real(DP) :: energy_i, distance

        overlap = .false.
        energy = 0._DP
        current => this%beginning%next
        if (.not. associated(current%next)) return
        do
            next => current%next
            if (.not. particle%same_type .or. particle%i /= current%i) then
                distance = this%periodic_box%distance(particle%position, current%position)
                call pair_potential%meet(overlap, energy_i, distance)
                if (overlap) return
                energy = energy + energy_i
            end if
            if (.not. associated(next%next)) return
            current => next
        end do
    end subroutine Abstract_Visitable_List_visit

    subroutine Abstract_Visitable_List_add(this, particle)
        class(Abstract_Visitable_List), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: particle

        type(Concrete_Linkable_Node), pointer :: previous => null(), new => null(), next => null()

        previous => this%beginning
        next => previous%next
        allocate(new)
        new%next => previous%next
        previous%next => new
        new%i = particle%i
        new%position = particle%position
    end subroutine Abstract_Visitable_List_add

    subroutine Abstract_Visitable_List_remove(this, i_particle)
        class(Abstract_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_particle

        type(Concrete_Linkable_Node), pointer :: previous => null(), current => null(), next => null()

        previous => this%beginning
        current => previous%next
        do
            next => current%next
            if (current%i == i_particle) then
                previous%next => current%next
                deallocate(current)
                current => next
                exit
            else
                previous => current
            end if
            current => next
        end do
    end subroutine Abstract_Visitable_List_remove

!end implementation Abstract_Visitable_List

!implementation Concrete_Visitable_Array

    subroutine Concrete_Visitable_Array_construct(this, periodic_box)
        class(Concrete_Visitable_Array), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box

        integer :: initial_size

        this%periodic_box => periodic_box
        this%num_nodes = 0
        initial_size = 10
        allocate(this%nodes(initial_size))
    end subroutine Concrete_Visitable_Array_construct

    subroutine Concrete_Visitable_Array_destroy(this)
        class(Concrete_Visitable_Array), intent(inout) :: this

        if (allocated(this%nodes)) deallocate(this%nodes)
    end subroutine Concrete_Visitable_Array_destroy

    subroutine Concrete_Visitable_Array_set(this, i_target, particle)
        class(Concrete_Visitable_Array), intent(inout) :: this
        integer, intent(in) :: i_target
        type(Concrete_Particle), intent(in) :: particle

        integer :: i_node

        do i_node = 1, this%num_nodes
            if (this%nodes(i_node)%i == i_target) then
                this%nodes(i_node)%i = particle%i
                this%nodes(i_node)%position = particle%position
                return
            end if
        end do
    end subroutine Concrete_Visitable_Array_set

    subroutine Concrete_Visitable_Array_visit(this, particle, pair_potential, overlap, energy)
        class(Concrete_Visitable_Array), intent(in) :: this
        type(Concrete_Particle), intent(in) :: particle
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy

        real(DP) :: energy_i, distance
        integer :: i_node

        overlap = .false.
        energy = 0._DP
        do i_node = 1, this%num_nodes
            if (.not. particle%same_type .or. particle%i /= this%nodes(i_node)%i) then
                distance = this%periodic_box%distance(particle%position, &
                    this%nodes(i_node)%position)
                call pair_potential%meet(overlap, energy_i, distance)
                if (overlap) return
                energy = energy + energy_i
            end if
        end do
    end subroutine Concrete_Visitable_Array_visit

    subroutine Concrete_Visitable_Array_add(this, particle)
        class(Concrete_Visitable_Array), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: particle

        this%num_nodes = this%num_nodes + 1
        if (size(this%nodes) < this%num_nodes) then
            call increase_nodes_size(this%nodes)
        end if
        this%nodes(this%num_nodes)%i = particle%i
        this%nodes(this%num_nodes)%position = particle%position
    end subroutine Concrete_Visitable_Array_add

    subroutine Concrete_Visitable_Array_remove(this, i_particle)
        class(Concrete_Visitable_Array), intent(inout) :: this
        integer, intent(in) :: i_particle

        type(Concrete_Particle) :: last

        if (i_particle /= this%nodes(this%num_nodes)%i) then
            last%i = this%nodes(this%num_nodes)%i
            last%position = this%nodes(this%num_nodes)%position
            call this%set(i_particle, last)
        end if
        this%num_nodes = this%num_nodes - 1
    end subroutine Concrete_Visitable_Array_remove

!end implementation Concrete_Visitable_Array

end module class_visitable_list
