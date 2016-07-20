module classes_box_size_checker

use procedures_errors, only: error_exit, warning_continue
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use classes_visitable_walls, only: Abstract_Visitable_Walls

implicit none

private

    type, abstract, public :: Abstract_Box_Size_Checker
    private
        class(Abstract_Reciprocal_Lattice), pointer :: reciprocal_lattice => null()
        class(Abstract_Visitable_Walls), pointer :: visitable_walls => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: check => Abstract_check
    end type Abstract_Box_Size_Checker

    type, extends(Abstract_Box_Size_Checker), public :: Concrete_Box_Size_Checker

    end type Concrete_Box_Size_Checker

    type, extends(Abstract_Box_Size_Checker), public :: Null_Box_Size_Checker
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: check => Null_check
    end type Null_Box_Size_Checker

contains

!implementation Abstract_Box_Size_Checker

    subroutine Abstract_construct(this, reciprocal_lattice, visitable_walls)
        class(Abstract_Box_Size_Checker), intent(out) :: this
        class(Abstract_Reciprocal_Lattice), target, intent(in) :: reciprocal_lattice
        class(Abstract_Visitable_Walls), target, intent(in) :: visitable_walls

        this%reciprocal_lattice => reciprocal_lattice
        this%visitable_walls => visitable_walls
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Box_Size_Checker), intent(out) :: this

        this%visitable_walls => null()
        this%reciprocal_lattice => null()
    end subroutine Abstract_destroy

    subroutine Abstract_check(this)
        class(Abstract_Box_Size_Checker), intent(in) :: this

        if (this%reciprocal_lattice%is_sparse_in_z()) then
            call warning_continue("Abstract_Box_Size_Checker: check: reciprocal_lattice is too "//&
                "sparse in z direction.")
        end if
        if (this%visitable_walls%are_outside_box()) then
            call error_exit("Abstract_Box_Size_Checker: check: "//&
                "visitable_walls are outside the box.")
        end if
    end subroutine Abstract_check

!end implementation Abstract_Box_Size_Checker

!implementation Null_Box_Size_Checker

    subroutine Null_construct(this, reciprocal_lattice, visitable_walls)
        class(Null_Box_Size_Checker), intent(out) :: this
        class(Abstract_Reciprocal_Lattice), target, intent(in) :: reciprocal_lattice
        class(Abstract_Visitable_Walls), target, intent(in) :: visitable_walls
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Box_Size_Checker), intent(out) :: this
    end subroutine Null_destroy

    subroutine Null_check(this)
        class(Null_Box_Size_Checker), intent(in) :: this
    end subroutine Null_check

!end implementation Null_Box_Size_Checker

end module classes_box_size_checker
