module procedures_boxes_size_checker_factory

use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use classes_visitable_walls, only: Abstract_Visitable_Walls
use classes_box_size_checker, only: Abstract_Box_Size_Checker, Concrete_Box_Size_Checker, &
    Null_Box_Size_Checker
use procedures_environment_inquirers, only: use_reciprocal_lattice, use_walls

implicit none

private
public :: create, destroy

contains

    !> @todo
    !> Replace all() by any()?
    subroutine create(box_size_checkers, reciprocal_lattices, visitable_walls)
        class(Abstract_Box_Size_Checker), allocatable, intent(out) :: box_size_checkers(:)
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattices(:)
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls(:)

        integer :: i_box

        if (all(use_reciprocal_lattice(reciprocal_lattices)) .or. all(use_walls(visitable_walls))) &
            then
            allocate(Concrete_Box_Size_Checker :: box_size_checkers(size(visitable_walls)))
        else
            allocate(Null_Box_Size_Checker :: box_size_checkers(size(visitable_walls)))
        end if
        do i_box = 1, size(box_size_checkers)
            call box_size_checkers(i_box)%construct(reciprocal_lattices(i_box), &
                visitable_walls(i_box))
        end do
    end subroutine create

    subroutine destroy(box_size_checkers)
        class(Abstract_Box_Size_Checker), allocatable, intent(inout) :: box_size_checkers(:)

        integer :: i_box

        if (allocated(box_size_checkers)) then
            do i_box = 1, size(box_size_checkers), -1
                call box_size_checkers(i_box)%destroy()
            end do
            deallocate(box_size_checkers)
        end if
    end subroutine destroy

end module procedures_boxes_size_checker_factory
