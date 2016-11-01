module procedures_line_writer_factory

use procedures_errors, only: error_exit
use types_string_wrapper, only: String_Wrapper
use classes_number_to_string, only: Concrete_Number_to_String
use procedures_environment_inquirers, only: box_size_can_change
use classes_external_field, only: Abstract_External_Field
use procedures_environment_inquirers, only: apply_external_field
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only: component_is_dipolar
use classes_pair_potential, only: Pair_Potential_Wrapper
use procedures_short_interactions_inquirers, only: component_interacts_with_wall
use classes_changed_box_size, only: Abstract_Changed_Box_Size
use classes_line_writer, only: Abstract_Line_Writer, Concrete_Line_Writer, Null_Line_Writer, &
    Line_Writer_Wrapper

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_teleportations_successes
    module procedure :: create_volumes_change_success
    module procedure :: create_walls_energies
    module procedure :: create_field_energies
    module procedure :: create_line
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_line
    module procedure :: destroy_teleportations_successes
end interface

contains

    !> @todo any(can_translate(:, i_box)): better condition?
    subroutine create_teleportations_successes(successes, make_directory_cmd, separator, directory,&
        can_translate)
        type(Line_Writer_Wrapper), allocatable, intent(out) :: successes(:, :)
        character(len=*), intent(in) :: make_directory_cmd, separator, directory
        logical, intent(in) :: can_translate(:, :)

        type(Concrete_Number_to_String) :: string
        integer :: num_boxes, i_box, j_box
        integer :: teleporation_stat

        num_boxes = size(can_translate, 2)
        if (num_boxes > 1 .and. any(can_translate)) then
            call execute_command_line(make_directory_cmd//" "//directory, &
                exitstat=teleporation_stat)
            if (teleporation_stat /= 0) call error_exit("procedures_line_writer_factory: "//&
                "create_teleportations_successes: "//directory//" directory can't be created.")
        end if

        allocate(successes(num_boxes, num_boxes))
        do j_box = 1, size(successes, 2)
            do i_box = 1, size(successes, 1)
                if (num_boxes > 1 .and. i_box /= j_box .and. any(can_translate(:, i_box))) then
                    allocate(Concrete_Line_Writer :: successes(i_box, j_box)%writer)
                else
                    allocate(Null_Line_Writer :: successes(i_box, j_box)%writer)
                end if
                call successes(i_box, j_box)%writer%construct(directory//separator//"successes_"//&
                    string%get(i_box)//"_to_"//string%get(j_box)//".out", can_translate(:, i_box))
            end do
        end do
    end subroutine create_teleportations_successes

    subroutine destroy_teleportations_successes(successes)
        type(Line_Writer_Wrapper), allocatable, intent(inout) :: successes(:, :)

        integer :: i_box, j_box

        if (allocated(successes)) then
            do j_box = size(successes, 2), 1, -1
                do i_box = size(successes, 1), 1, -1
                    call successes(i_box, j_box)%writer%destroy()
                end do
            end do
            deallocate(successes)
        end if
    end subroutine destroy_teleportations_successes

    subroutine create_volumes_change_success(successes, filename, changed_boxes_size)
        class(Abstract_Line_Writer), allocatable, intent(out) :: successes
        character(len=*), intent(in) :: filename
        class(Abstract_Changed_Box_Size), intent(in) :: changed_boxes_size(:)

        logical :: selector(size(changed_boxes_size))

        selector = box_size_can_change(changed_boxes_size)

        if (any(selector)) then
            allocate(Concrete_Line_Writer :: successes)
        else
            allocate(Null_Line_Writer :: successes)
        end if
        call successes%construct(filename, selector)
    end subroutine create_volumes_change_success

    subroutine create_walls_energies(energies, filename, wall_pairs, visit_energies)
        class(Abstract_Line_Writer), allocatable, intent(out) :: energies
        character(len=*), intent(in) :: filename
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        logical, intent(in) :: visit_energies

        logical :: selector(size(wall_pairs))
        integer :: i_component

        do i_component = 1, size(wall_pairs)
            selector(i_component) = component_interacts_with_wall(wall_pairs(i_component)%potential)
        end do

        if (any(selector) .and. visit_energies) then
            allocate(Concrete_Line_Writer :: energies)
        else
            allocate(Null_Line_Writer :: energies)
        end if
        call energies%construct(filename, selector)
    end subroutine create_walls_energies

    subroutine create_field_energies(energies, filename, external_field, components, visit_energies)
        class(Abstract_Line_Writer), allocatable, intent(out) :: energies
        character(len=*), intent(in) :: filename
        class(Abstract_External_Field), intent(in) :: external_field
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: visit_energies

        logical :: selector(size(components))
        integer :: i_component

        do i_component = 1, size(selector)
            selector(i_component) = component_is_dipolar(components(i_component)%dipole_moments)
        end do

        if (any(selector) .and. apply_external_field(external_field) .and. visit_energies) then
            allocate(Concrete_Line_Writer :: energies)
        else
            allocate(Null_Line_Writer :: energies)
        end if
        call energies%construct(filename, selector)
    end subroutine create_field_energies

    !> @note Condition may not be optimal.
    subroutine create_line(writers, paths, filename, selectors)
        class(Abstract_Line_Writer), allocatable, intent(out) :: writers(:)
        type(String_Wrapper), intent(in) :: paths(:)
        character(len=*), intent(in) :: filename
        logical, intent(in) :: selectors(:, :)

        integer :: i_element

        if (any(selectors)) then
            allocate(Concrete_Line_Writer :: writers(size(paths)))
        else
            allocate(Null_Line_Writer :: writers(size(paths)))
        end if

        do i_element = 1, size(writers)
            call writers(i_element)%construct(paths(i_element)%string//filename, &
                selectors(:, i_element))
        end do
    end subroutine create_line

    subroutine destroy_line(writers)
        class(Abstract_Line_Writer), allocatable, intent(inout) :: writers(:)

        integer :: i_element

        if (allocated(writers)) then
            do i_element = size(writers), 1, -1
                call writers(i_element)%destroy()
            end do
            deallocate(writers)
        end if
    end subroutine destroy_line

    subroutine destroy_element(writer)
        class(Abstract_Line_Writer), allocatable, intent(inout) :: writer

        if (allocated(writer)) then
            call writer%destroy()
            deallocate(writer)
        end if
    end subroutine destroy_element

end module procedures_line_writer_factory

