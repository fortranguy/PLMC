module procedures_cells_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use types_component_wrapper, only: Component_Wrapper
use classes_hard_contact, only: Abstract_Hard_Contact
use classes_pair_potential, only: Pair_Potential_Line
use procedures_neighbour_cells_factory, only: neighbour_cells_create_triangle => create_triangle, &
    neighbour_cells_create_element => create_element, neighbour_cells_destroy_triangle => &
    destroy_triangle,  neighbour_cells_destroy_element => destroy_element, &
    neighbour_cells_allocate_triangle => allocate_triangle
use classes_visitable_list, only: Abstract_Visitable_List
use procedures_visitable_cells_factory, only: visitable_cells_create => create, &
    visitable_cells_destroy => destroy
use procedures_visitable_cells_memento_factory, only: visitable_cells_memento_create => create, &
    visitable_cells_memento_destroy => destroy
use types_cells_wrapper, only: Cells_Wrapper

implicit none

private
public :: create, destroy, allocate_triangle

interface create
    module procedure :: create_wrappers
    module procedure :: visitable_cells_create
    module procedure :: visitable_cells_memento_create
    module procedure :: neighbour_cells_create_triangle, neighbour_cells_create_element
end interface create

interface destroy
    module procedure :: neighbour_cells_destroy_triangle, neighbour_cells_destroy_element
    module procedure :: visitable_cells_memento_destroy
    module procedure :: visitable_cells_destroy
    module procedure :: destroy_wrappers
end interface destroy

interface allocate_triangle
    module procedure :: neighbour_cells_allocate_triangle
end interface allocate_triangle

contains

    subroutine create_wrappers(cells, periodic_boxes, accessible_domains, components, hard_contact,&
        components_pairs, list_mold, interact)
        type(Cells_Wrapper), allocatable, intent(out) :: cells(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        class(Abstract_Parallelepiped_Domain), intent(in) :: accessible_domains(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        class(Abstract_Hard_Contact), intent(in) :: hard_contact
        type(Pair_Potential_Line), intent(in) :: components_pairs(:)
        class(Abstract_Visitable_List), intent(in) :: list_mold
        logical, intent(in) :: interact

        integer :: i_box

        allocate(cells(size(periodic_boxes)))
        do i_box = 1, size(cells)
            call neighbour_cells_create_triangle(cells(i_box)%neighbour_cells, &
                periodic_boxes(i_box), accessible_domains(i_box), hard_contact, components_pairs, &
                interact)
            call visitable_cells_create(cells(i_box)%visitable_cells, periodic_boxes(i_box), &
                components(:, i_box), hard_contact, components_pairs, cells(i_box)%neighbour_cells,&
                list_mold, interact)
        end do
    end subroutine create_wrappers

    subroutine destroy_wrappers(cells)
        type(Cells_Wrapper), allocatable, intent(inout) :: cells(:)

        integer :: i_box

        if (allocated(cells)) then
            do i_box = size(cells), 1, -1
                call visitable_cells_destroy(cells(i_box)%visitable_cells)
                call neighbour_cells_destroy_triangle(cells(i_box)%neighbour_cells)
            end do
            deallocate(cells)
        end if
    end subroutine destroy_wrappers


end module procedures_cells_factory
