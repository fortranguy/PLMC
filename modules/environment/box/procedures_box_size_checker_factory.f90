module procedures_box_size_checker_factory

use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use classes_visitable_walls, only: Abstract_Visitable_Walls
use classes_box_size_checker, only: Abstract_Box_Size_Checker, Concrete_Box_Size_Checker, &
    Null_Box_Size_Checker
use procedures_environment_inquirers, only: use_reciprocal_lattice, use_walls

implicit none

private
public :: create, destroy

contains

    subroutine create(box_size_checker, reciprocal_lattice, visitable_walls)
        class(Abstract_Box_Size_Checker), allocatable, intent(out) :: box_size_checker
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls

        if (use_reciprocal_lattice(reciprocal_lattice) .or. use_walls(visitable_walls)) then
            allocate(Concrete_Box_Size_Checker :: box_size_checker)
        else
            allocate(Null_Box_Size_Checker :: box_size_checker)
        end if
        call box_size_checker%construct(reciprocal_lattice, visitable_walls)
    end subroutine create

    subroutine destroy(box_size_checker)
        class(Abstract_Box_Size_Checker), allocatable, intent(inout) :: box_size_checker

        if (allocated(box_size_checker)) then
            call box_size_checker%destroy()
            deallocate(box_size_checker)
        end if
    end subroutine destroy

end module procedures_box_size_checker_factory
