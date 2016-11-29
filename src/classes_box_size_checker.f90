module classes_box_size_checker

use procedures_errors, only: error_exit, warning_continue
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use classes_visitable_walls, only: Abstract_Visitable_Walls

implicit none

private

    type, abstract, public :: Abstract_Box_Size_Checker
    private
        class(Abstract_Parallelepiped_Domain), pointer :: accessible_domain => null()
        class(Abstract_Parallelepiped_Domain), pointer :: field_domain => null()
        class(Abstract_Reciprocal_Lattice), pointer :: reciprocal_lattice => null()
        class(Abstract_Visitable_Walls), pointer :: visitable_walls => null()
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_check), deferred :: check
        procedure, private :: check_core => Abstract_check_core
    end type Abstract_Box_Size_Checker

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Box_Size_Checker
            class(Abstract_Box_Size_Checker), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_check(this)
        import :: Abstract_Box_Size_Checker
            class(Abstract_Box_Size_Checker), intent(in) :: this
        end subroutine Abstract_check

    end interface

    type, extends(Abstract_Box_Size_Checker), public :: Generating_Box_Size_Checker
    contains
        procedure :: construct => Generating_construct
        procedure :: destroy => Generating_destroy
        procedure :: check => Generating_check
    end type Generating_Box_Size_Checker

    type, extends(Abstract_Box_Size_Checker), public :: Exploring_Box_Size_Checker
    private
        class(Abstract_Parallelepiped_Domain), pointer :: particle_insertion_domain => null()
    contains
        procedure :: construct => Exploring_construct
        procedure :: destroy => Exploring_destroy
        procedure :: check => Exploring_check
    end type Exploring_Box_Size_Checker

contains

!implementation Abstract_Box_Size_Checker

    subroutine Abstract_check_core(this)
        class(Abstract_Box_Size_Checker), intent(in) :: this

        if (.not.this%accessible_domain%is_boxed()) then
            call error_exit("Abstract_Box_Size_Checker: accessible_domain is not boxed.")
        end if
        if (.not.this%field_domain%is_boxed()) then
            call error_exit("Abstract_Box_Size_Checker: field_domain is not boxed.")
        end if
        if (this%reciprocal_lattice%is_sparse_in_z()) then
            call warning_continue("Abstract_Box_Size_Checker: check: reciprocal_lattice is too "//&
                "sparse in z direction.")
        end if
        if (this%visitable_walls%are_outside_box()) then
            call error_exit("Abstract_Box_Size_Checker: check: "//&
                "visitable_walls are outside the box.")
        end if
    end subroutine Abstract_check_core

!end implementation Abstract_Box_Size_Checker

!implementation Generating_Box_Size_Checker

    subroutine Generating_construct(this, accessible_domain, field_domain, reciprocal_lattice, &
        visitable_walls)
        class(Generating_Box_Size_Checker), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: accessible_domain
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: field_domain
        class(Abstract_Reciprocal_Lattice), target, intent(in) :: reciprocal_lattice
        class(Abstract_Visitable_Walls), target, intent(in) :: visitable_walls

        this%accessible_domain => accessible_domain
        this%field_domain => field_domain
        this%reciprocal_lattice => reciprocal_lattice
        this%visitable_walls => visitable_walls
    end subroutine Generating_construct

    subroutine Generating_destroy(this)
        class(Generating_Box_Size_Checker), intent(inout) :: this

        this%visitable_walls => null()
        this%reciprocal_lattice => null()
        this%field_domain => null()
        this%accessible_domain => null()
    end subroutine Generating_destroy

    subroutine Generating_check(this)
        class(Generating_Box_Size_Checker), intent(in) :: this

        call this%check_core()
    end subroutine Generating_check

!end implementation Generating_Box_Size_Checker

!implementation Exploring_Box_Size_Checker

    subroutine Exploring_construct(this, accessible_domain, field_domain, reciprocal_lattice, &
        visitable_walls, particle_insertion_domain)
        class(Exploring_Box_Size_Checker), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: accessible_domain
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: field_domain
        class(Abstract_Reciprocal_Lattice), target, intent(in) :: reciprocal_lattice
        class(Abstract_Visitable_Walls), target, intent(in) :: visitable_walls
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: particle_insertion_domain

        this%accessible_domain => accessible_domain
        this%field_domain => field_domain
        this%reciprocal_lattice => reciprocal_lattice
        this%visitable_walls => visitable_walls
        this%particle_insertion_domain => particle_insertion_domain
    end subroutine Exploring_construct

    subroutine Exploring_destroy(this)
        class(Exploring_Box_Size_Checker), intent(inout) :: this

        this%particle_insertion_domain => null()
        this%visitable_walls => null()
        this%reciprocal_lattice => null()
        this%field_domain => null()
        this%accessible_domain => null()
    end subroutine Exploring_destroy

    subroutine Exploring_check(this)
        class(Exploring_Box_Size_Checker), intent(in) :: this

        call this%check_core()
        if (.not.this%particle_insertion_domain%is_boxed()) then
            call error_exit("Exploring_Box_Size_Checker: particle_insertion_domain is not boxed.")
        end if
    end subroutine Exploring_check

!end implementation Exploring_Box_Size_Checker

end module classes_box_size_checker
