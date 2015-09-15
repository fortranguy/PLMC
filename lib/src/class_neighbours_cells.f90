module class_neighbour_cells

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_visitable_list, only: Abstract_Visitable_List
use module_particles, only: Concrete_Particle

implicit none

private

    type, abstract, public :: Abstract_Neighbour_Cells
    private
        class(Abstract_Visitable_List), pointer :: linked_list
    contains
        procedure :: construct => Abstract_Neighbour_Cells_construct
        procedure :: destroy => Abstract_Neighbour_Cells_destroy
        procedure :: visit => Abstract_Neighbour_Cells_visit
    end type Abstract_Neighbour_Cells

    type, extends(Abstract_Neighbour_Cells), public :: Concrete_Neighbour_Cells

    end type Concrete_Neighbour_Cells

contains

    subroutine Abstract_Neighbour_Cells_construct(this, linked_list)
        class(Abstract_Neighbour_Cells), intent(out) :: this
        class(Abstract_Visitable_List), target, intent(in) :: linked_list

        this%linked_list => linked_list
    end subroutine Abstract_Neighbour_Cells_construct

    subroutine Abstract_Neighbour_Cells_destroy(this)
        class(Abstract_Neighbour_Cells), intent(inout) :: this

        this%linked_list => null()
    end subroutine Abstract_Neighbour_Cells_destroy

    subroutine Abstract_Neighbour_Cells_visit(this, particle, energy)
        class(Abstract_Neighbour_Cells), intent(in) :: this
        type(Concrete_Particle), intent(in) :: particle
        real(DP), intent(out) :: energy

        !call this%linked_list%visit(particle, energy)
    end subroutine Abstract_Neighbour_Cells_visit

end module class_neighbour_cells
