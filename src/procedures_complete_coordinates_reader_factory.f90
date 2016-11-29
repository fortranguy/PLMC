module procedures_complete_coordinates_reader_factory

use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use procedures_boxes_size_checker_factory, only: boxes_size_checker_create => create, &
    boxes_size_checker_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use classes_box_size_checker, only: Abstract_Box_Size_Checker
use types_component_wrapper, only: Component_Wrapper
use classes_component_coordinates_reader, only: Component_Coordinates_Reader_wrapper
use procedures_component_coordinates_reader_factory, only: &
    component_coordinates_reader_create => create, component_coordinates_reader_destroy => destroy
use classes_complete_coordinates_reader, only: Abstract_Complete_Coordinates_Reader, &
    Concrete_Complete_Coordinates_Reader

implicit none

private
public :: create, destroy

contains

    subroutine create(coordinates, environment, components, particle_insertion_domains)
        class(Abstract_Complete_Coordinates_Reader), allocatable, intent(out) :: coordinates
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:, :)
        class(Abstract_Parallelepiped_Domain), optional, intent(in) ::particle_insertion_domains(:)

        class(Abstract_Box_Size_Checker), allocatable :: boxes_size_checker(:)
        type(Component_Coordinates_Reader_wrapper), allocatable :: components_coordinates(:, :)

        if (present(particle_insertion_domains)) then
            call boxes_size_checker_create(boxes_size_checker, environment%accessible_domains, &
                environment%fields_domain, environment%reciprocal_lattices, environment%&
                visitable_walls, particle_insertion_domains)
        else
            call boxes_size_checker_create(boxes_size_checker, environment%accessible_domains, &
                environment%fields_domain, environment%reciprocal_lattices, environment%&
                visitable_walls)
        end if
        call component_coordinates_reader_create(components_coordinates, components)

        allocate(Concrete_Complete_Coordinates_Reader :: coordinates)
        call coordinates%construct(environment%periodic_boxes, boxes_size_checker, &
            components_coordinates)

        call component_coordinates_reader_destroy(components_coordinates)
        call boxes_size_checker_destroy(boxes_size_checker)
    end subroutine create

    subroutine destroy(coordinates)
        class(Abstract_Complete_Coordinates_Reader), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy

end module procedures_complete_coordinates_reader_factory
