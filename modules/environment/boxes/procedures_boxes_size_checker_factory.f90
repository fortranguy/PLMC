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
    subroutine create(boxes_size_checker, reciprocal_lattices, visitable_walls)
        class(Abstract_Box_Size_Checker), allocatable, intent(out) :: boxes_size_checker(:)
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattices(:)
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls(:)

        integer :: i_box

        if (all(use_reciprocal_lattice(reciprocal_lattices)) .or. all(use_walls(visitable_walls))) &
            then
            allocate(Concrete_Box_Size_Checker :: boxes_size_checker(size(visitable_walls)))
        else
            allocate(Null_Box_Size_Checker :: boxes_size_checker(size(visitable_walls)))
        end if
        do i_box = 1, size(boxes_size_checker)
            call boxes_size_checker(i_box)%construct(reciprocal_lattices(i_box), &
                visitable_walls(i_box))
        end do
    end subroutine create

    subroutine destroy(boxes_size_checker)
        class(Abstract_Box_Size_Checker), allocatable, intent(inout) :: boxes_size_checker(:)

        integer :: i_box

        if (allocated(boxes_size_checker)) then
            do i_box = size(boxes_size_checker), 1, -1
                call boxes_size_checker(i_box)%destroy()
            end do
            deallocate(boxes_size_checker)
        end if
    end subroutine destroy

end module procedures_boxes_size_checker_factory
