module procedures_boxes_size_checker_factory

use procedures_errors, only: error_exit
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use classes_visitable_walls, only: Abstract_Visitable_Walls
use classes_box_size_checker, only: Abstract_Box_Size_Checker, Generating_Box_Size_Checker, &
    Exploring_Box_Size_Checker

implicit none

private
public :: create, destroy

contains

    subroutine create(boxes_size_checker, accessible_domains, fields_domain, reciprocal_lattices, &
        visitable_walls, particle_insertion_domains)
        class(Abstract_Box_Size_Checker), allocatable, intent(out) :: boxes_size_checker(:)
        class(Abstract_Parallelepiped_Domain), intent(in) :: accessible_domains(:)
        class(Abstract_Parallelepiped_Domain), intent(in) :: fields_domain(:)
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattices(:)
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls(:)
        class(Abstract_Parallelepiped_Domain), optional, intent(in) :: particle_insertion_domains(:)

        integer :: i_box

        if (present(particle_insertion_domains)) then
            allocate(Exploring_Box_Size_Checker :: boxes_size_checker(size(accessible_domains)))
        else
            allocate(Generating_Box_Size_Checker :: boxes_size_checker(size(accessible_domains)))
        end if

        select type (boxes_size_checker)
            type is (Generating_Box_Size_Checker)
                do i_box = 1, size(boxes_size_checker)
                    call boxes_size_checker(i_box)%construct(accessible_domains(i_box), &
                        fields_domain(i_box), reciprocal_lattices(i_box), visitable_walls(i_box))
                end do
            type is (Exploring_Box_Size_Checker)
                do i_box = 1, size(boxes_size_checker)
                    call boxes_size_checker(i_box)%construct(accessible_domains(i_box), &
                        fields_domain(i_box), reciprocal_lattices(i_box), visitable_walls(i_box), &
                        particle_insertion_domains(i_box))
                end do
            class default
                call error_exit("procedures_boxes_size_checker_factory: create: "//&
                    "boxes_size_checker: type unknown.")
        end select
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
