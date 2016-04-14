module class_box_size_checker

use procedures_errors, only: error_exit, warning_continue
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use class_walls_potential, only: Abstract_Walls_Potential

implicit none

private

    type, abstract, public :: Abstract_Box_Size_Checker
    private
        class(Abstract_Reciprocal_Lattice), pointer :: reciprocal_lattice => null()
        class(Abstract_Walls_Potential), pointer :: walls_potential => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: check => Abstract_check
    end type Abstract_Box_Size_Checker

contains

    subroutine Abstract_construct(this, reciprocal_lattice, walls_potential)
        class(Abstract_Box_Size_Checker), intent(out) :: this
        class(Abstract_Reciprocal_Lattice), target, intent(in) :: reciprocal_lattice
        class(Abstract_Walls_Potential), target, intent(in) :: walls_potential

        this%reciprocal_lattice => reciprocal_lattice
        this%walls_potential => walls_potential
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Box_Size_Checker), intent(out) :: this

        this%walls_potential => null()
        this%reciprocal_lattice => null()
    end subroutine Abstract_destroy

    subroutine Abstract_check(this)
        class(Abstract_Box_Size_Checker), intent(in) :: this

        if (this%reciprocal_lattice%is_sparse_in_z()) then
            call warning_continue("Abstract_Box_Size_Checker: check: reciprocal_lattice is too "//&
                "sparse in z direction.")
        end if
        if (this%walls_potential%are_outside_box()) then
            call error_exit("Abstract_Box_Size_Checker: check: walls are outside the box.")
        end if
    end subroutine Abstract_check

end module class_box_size_checker
