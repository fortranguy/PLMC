module class_visitable_list

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use module_particles, only: Concrete_Particle
use module_nodes, only: Concrete_Node, &
    deallocate_list
use class_periodic_box, only: Abstract_Periodic_Box
use class_positions, only: Abstract_Positions
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Visitable_List
    private
        type(Concrete_Node), pointer :: beginning
        class(Abstract_Periodic_Box), pointer :: periodic_box
        class(Abstract_Positions), pointer :: positions
        class(Abstract_Pair_Potential), pointer :: pair_potential
    contains
        procedure :: construct => Abstract_Visitable_List_construct
        procedure :: destroy => Abstract_Visitable_List_destroy
        procedure :: set => Abstract_Visitable_List_set
        procedure :: visit => Abstract_Visitable_List_visit
        procedure :: add => Abstract_Visitable_List_add
        procedure :: remove => Abstract_Visitable_List_remove
        procedure :: overwrite => Abstract_Visitable_List_overwrite
    end type Abstract_Visitable_List

    type, extends(Abstract_Visitable_List), public :: Concrete_Visitable_List

    end type Concrete_Visitable_List

contains

    subroutine Abstract_Visitable_List_construct(this, periodic_box, positions)
        class(Abstract_Visitable_List), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Positions), target, intent(in) :: positions

        type(Concrete_Node), pointer :: current => null(), next => null()

        this%periodic_box => periodic_box
        this%positions => positions

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
        this%pair_potential => null()
        this%positions => null()
        this%periodic_box => null()
    end subroutine Abstract_Visitable_List_destroy

    subroutine Abstract_Visitable_List_set(this, pair_potential)
        class(Abstract_Visitable_List), intent(inout) :: this
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential

        this%pair_potential => pair_potential
    end subroutine Abstract_Visitable_List_set

    subroutine Abstract_Visitable_List_visit(this, same_type, particle, overlap, energy)
        class(Abstract_Visitable_List), intent(in) :: this
        logical, intent(in) :: same_type
        type(Concrete_Particle), intent(in) :: particle
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy

        type(Concrete_Node), pointer :: current => null(), next => null()
        real(DP) :: energy_i, distance

        overlap = .false.
        energy = 0._DP
        current => this%beginning%next
        do
            next => current%next
            if (.not. same_type .or. particle%i /= current%i) then
                distance = this%periodic_box%distance(particle%position, &
                                                      this%positions%get(current%i))
                call this%pair_potential%meet(distance, overlap, energy_i)
                if (overlap) return
                energy = energy + energy_i
            end if
            if (.not. associated(next%next)) exit
            current => next
        end do
    end subroutine Abstract_Visitable_List_visit

    subroutine Abstract_Visitable_List_add(this, i_particle)
        class(Abstract_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_particle

        type(Concrete_Node), pointer :: previous => null(), new => null(), next => null()

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

        type(Concrete_Node), pointer :: previous => null(), current => null(), next => null()

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

    subroutine Abstract_Visitable_List_overwrite(this, i_target, i_value)
        class(Abstract_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_target, i_value

        type(Concrete_Node), pointer :: previous => null(), current => null(), next => null()

        previous => this%beginning
        current => previous%next
        do
            next => current%next
            if (current%i == i_target) then
                current%i = i_value
                exit
            else
                previous => current
            end if
            current => next
        end do
    end subroutine Abstract_Visitable_List_overwrite

end module class_visitable_list
