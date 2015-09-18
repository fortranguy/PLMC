module class_cells_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_pair_potential, only: Abstract_Pair_Potential
use module_particles, only: Concrete_Particle
use class_visitable_cells, only: Abstract_Visitable_Cells

implicit none

private

    type, public :: Cells_Potential_Facade
    private
        class(Abstract_Pair_Potential), pointer :: pair_potential
        class(Abstract_Visitable_Cells), pointer :: visitable_cells
    contains
        procedure :: construct => Cells_Potential_Facade_construct
        procedure :: destroy => Cells_Potential_Facade_destroy
        procedure :: set => Cells_Potential_Facade_set
        procedure :: visit => Cells_Potential_Facade_visit
    end type Cells_Potential_Facade

contains

    subroutine Cells_Potential_Facade_construct(this, visitable_cells)
        class(Cells_Potential_Facade), intent(out) :: this
        class(Abstract_Visitable_Cells), target, intent(in) :: visitable_cells

        this%visitable_cells => visitable_cells
    end subroutine Cells_Potential_Facade_construct

    subroutine Cells_Potential_Facade_destroy(this)
        class(Cells_Potential_Facade), intent(inout) :: this

        this%visitable_cells => null()
        this%pair_potential => null()
    end subroutine Cells_Potential_Facade_destroy

    subroutine Cells_Potential_Facade_set(this, pair_potential)
        class(Cells_Potential_Facade), intent(inout) :: this
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential

        this%pair_potential => pair_potential
    end subroutine Cells_Potential_Facade_set

    subroutine Cells_Potential_Facade_visit(this, particle, overlap, energy)
        class(Cells_Potential_Facade), intent(in) :: this
        type(Concrete_Particle), intent(in) :: particle
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy

        call this%visitable_cells%visit(particle, this%pair_potential, overlap, energy)
    end subroutine Cells_Potential_Facade_visit

end module class_cells_potential
