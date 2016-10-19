module procedures_component_coordinates_reader_factory

use procedures_errors, only: error_exit
use classes_num_particles, only: Abstract_Num_Particles
use classes_component_coordinates, only: Abstract_Component_Coordinates
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only: component_has_positions, component_has_orientations
use types_component_coordinates_reader_selector, only: Component_Coordinates_Reader_Selector
use classes_component_coordinates_reader, only: Abstract_Component_Coordinates_Reader, &
    Concrete_Component_Coordinates_Reader, Concrete_Component_Positions_Reader, &
    Concrete_Component_Orientations_Reader, Null_Component_Coordinates_Reader, &
    Component_Coordinates_Reader_wrapper

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_rectangle
    module procedure :: create_element
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_rectangle
end interface destroy

contains

    subroutine create_rectangle(components_coordinates, components)
        type(Component_Coordinates_Reader_wrapper), allocatable, intent(out) :: &
            components_coordinates(:, :)
        type(Component_Wrapper), intent(in) :: components(:, :)

        type(Component_Coordinates_Reader_Selector) :: selector_i
        integer :: i_box, i_component

        allocate(components_coordinates(size(components, 1), size(components, 2)))
        do i_box = 1, size(components_coordinates, 2)
            do i_component = 1, size(components_coordinates, 1)
                associate(num_particles_i => components(i_component, i_box)%num_particles, &
                    positions_i => components(i_component, i_box)%positions, &
                    orientations_i => components(i_component, i_box)%orientations)
                    selector_i%read_positions = component_has_positions(positions_i)
                    selector_i%read_orientations = component_has_orientations(orientations_i)
                    call create(components_coordinates(i_component, i_box)%reader, num_particles_i,&
                        positions_i, orientations_i, selector_i)
                end associate
            end do
        end do
    end subroutine create_rectangle

    subroutine destroy_rectangle(components_coordinates)
        type(Component_Coordinates_Reader_wrapper), allocatable, intent(inout) :: &
            components_coordinates(:, :)

        integer :: i_box, i_component

        if (allocated(components_coordinates)) then
            do i_box = size(components_coordinates, 2), 1, -1
                do i_component = size(components_coordinates, 1), 1, -1
                    call destroy(components_coordinates(i_component, i_box)%reader)
                end do
            end do
            deallocate(components_coordinates)
        end if
    end subroutine destroy_rectangle

    subroutine create_element(coordinates, num_particles, positions, orientations, selector)
        class(Abstract_Component_Coordinates_Reader), allocatable, intent(out) :: coordinates
        class(Abstract_Num_Particles), intent(in) :: num_particles
        class(Abstract_Component_Coordinates), intent(in) :: positions, orientations
        type(Component_Coordinates_Reader_Selector), intent(in) :: selector

        call allocate_element(coordinates, selector)
        call construct_element(coordinates, num_particles, positions, orientations)
    end subroutine create_element

    subroutine allocate_element(coordinates, selector)
        class(Abstract_Component_Coordinates_Reader), allocatable, intent(out) :: coordinates
        type(Component_Coordinates_Reader_Selector), intent(in) :: selector

        if (selector%read_positions .and. selector%read_orientations) then
            allocate(Concrete_Component_Coordinates_Reader :: coordinates)
        else if (selector%read_positions .and. .not.selector%read_orientations) then
            allocate(Concrete_Component_Positions_Reader :: coordinates)
        else if (.not.selector%read_positions .and. selector%read_orientations) then
            allocate(Concrete_Component_Orientations_Reader :: coordinates)
        else
            allocate(Null_Component_Coordinates_Reader :: coordinates)
        end if
    end subroutine allocate_element

    subroutine construct_element(coordinates, num_particles, positions, orientations)
        class(Abstract_Component_Coordinates_Reader), intent(inout) :: coordinates
        class(Abstract_Num_Particles), intent(in) :: num_particles
        class(Abstract_Component_Coordinates), intent(in) :: positions, orientations

        select type (coordinates)
            type is (Concrete_Component_Coordinates_Reader)
                call coordinates%construct(num_particles, positions, orientations)
            type is (Concrete_Component_Positions_Reader)
                call coordinates%construct(num_particles, positions)
            type is (Concrete_Component_Orientations_Reader)
                call coordinates%construct(num_particles, orientations)
            type is (Null_Component_Coordinates_Reader)
                call coordinates%construct()
            class default
                call error_exit("procedures_component_coordinates_reader_factory: construct: "//&
                    "coordinates type unknown.")
        end select
    end subroutine construct_element

    subroutine destroy_element(coordinates)
        class(Abstract_Component_Coordinates_Reader), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy_element

end module procedures_component_coordinates_reader_factory
