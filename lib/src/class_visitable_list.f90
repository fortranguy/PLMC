module class_visitable_list

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use module_particles, only: Concrete_Particle
use module_nodes, only: Concrete_Node, &
    deallocate_list
use class_periodic_box, only: Abstract_Periodic_Box
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Visitable_List
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box
        type(Concrete_Node), pointer :: beginning
    contains
        procedure :: construct => Abstract_Visitable_List_construct
        procedure :: destroy => Abstract_Visitable_List_destroy
        procedure :: visit => Abstract_Visitable_List_visit
        procedure :: allocate => Abstract_Visitable_List_allocate
        procedure :: deallocate => Abstract_Visitable_List_deallocate
        procedure :: overwrite => Abstract_Visitable_List_overwrite
    end type Abstract_Visitable_List

    type, extends(Abstract_Visitable_List), public :: Concrete_Visitable_List

    end type Concrete_Visitable_List

contains

    subroutine Abstract_Visitable_List_construct(this, periodic_box)
        class(Abstract_Visitable_List), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box

        type(Concrete_Node), pointer :: current => null(), next => null()

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

    subroutine Abstract_Visitable_List_visit(this, particle, pair_potential, overlap, energy)
        class(Abstract_Visitable_List), intent(in) :: this
        type(Concrete_Particle), intent(in) :: particle
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy

        type(Concrete_Node), pointer :: current => null(), next => null()
        real(DP) :: energy_i, distance

        overlap = .false.
        energy = 0._DP
        current => this%beginning%next
        if (.not. associated(current%next)) return
        do
            next => current%next
            if (.not. particle%same_type .or. particle%i /= current%i) then
                distance = this%periodic_box%distance(particle%position, current%position)
                call pair_potential%meet(distance, overlap, energy_i)
                if (overlap) return
                energy = energy + energy_i
            end if
            if (.not. associated(next%next)) return
            current => next
        end do
    end subroutine Abstract_Visitable_List_visit

    subroutine Abstract_Visitable_List_allocate(this, particle)
        class(Abstract_Visitable_List), intent(inout) :: this
        type(Concrete_Particle), intent(in) :: particle

        type(Concrete_Node), pointer :: previous => null(), new => null(), next => null()

        previous => this%beginning
        next => previous%next
        allocate(new)
        new%next => previous%next
        previous%next => new
        new%i = particle%i
        new%position = particle%position
    end subroutine Abstract_Visitable_List_allocate

    subroutine Abstract_Visitable_List_deallocate(this, i_particle)
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
    end subroutine Abstract_Visitable_List_deallocate

    subroutine Abstract_Visitable_List_overwrite(this, i_target, value)
        class(Abstract_Visitable_List), intent(inout) :: this
        integer, intent(in) :: i_target
        type(Concrete_Particle), intent(in) :: value

        type(Concrete_Node), pointer :: previous => null(), current => null(), next => null()

        previous => this%beginning
        current => previous%next
        do
            next => current%next
            if (current%i == i_target) then
                current%i = value%i
                current%position = value%position
                exit
            else
                previous => current
            end if
            current => next
        end do
    end subroutine Abstract_Visitable_List_overwrite

end module class_visitable_list
