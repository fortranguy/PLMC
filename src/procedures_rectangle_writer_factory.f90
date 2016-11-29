module procedures_rectangle_writer_factory

use procedures_errors, only: error_exit
use types_string_wrapper, only: String_Wrapper
use classes_number_to_string, only: Concrete_Number_to_String
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only: component_can_exchange
use classes_rectangle_writer, only: Abstract_Rectangle_Writer, Concrete_Rectangle_Writer, &
    Null_Rectangle_Writer, Rectangle_Writer_Wrapper

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_swaps_successes
    module procedure :: create_transmutations_successes
end interface create

interface destroy
    module procedure :: destroy_line
    module procedure :: destroy_rectangle
end interface destroy

contains

    subroutine create_swaps_successes(successes, make_directory_cmd, separator, directory, &
        can_translate)
        type(Rectangle_Writer_Wrapper), allocatable, intent(out) :: successes(:, :)
        character(len=*), intent(in) :: make_directory_cmd, separator, directory
        logical, intent(in) :: can_translate(:, :)

        type(Concrete_Number_to_String) :: string
        integer :: swap_stat
        integer :: num_boxes, i_box, j_box
        logical :: selectors(size(can_translate, 1), size(can_translate, 1))
        integer :: i_component

        num_boxes = size(can_translate, 2)
        if (num_boxes > 1 .and. any(can_translate)) then
            swap_stat = 1
            call execute_command_line(make_directory_cmd//" "//directory, exitstat=swap_stat)
            if (swap_stat /= 0) call error_exit("procedures_rectangle_writer_factory: "//&
                "create_swaps_successes: "//directory//" directory can't be created.")
        end if

        selectors = .true.
        forall (i_component=1:size(can_translate, 1)) selectors(i_component, i_component) = .false.

        allocate(successes(num_boxes, num_boxes))
        do j_box = 1, size(successes, 2)
            do i_box = 1, size(successes, 1)
                if (i_box /= j_box .and. any(can_translate(:, i_box))) then
                    allocate(Concrete_Rectangle_Writer :: successes(i_box, j_box)%writer)
                else
                    allocate(Null_Rectangle_Writer :: successes(i_box, j_box)%writer)
                end if
                call successes(i_box, j_box)%writer%construct(directory//separator//"successes_"//&
                    string%get(i_box)//"_to_"//string%get(j_box)//".out", selectors)
            end do
        end do
    end subroutine create_swaps_successes

    subroutine destroy_rectangle(rectangles)
        type(Rectangle_Writer_Wrapper), allocatable, intent(inout) :: rectangles(:, :)

        integer :: i_element, j_element

        if (allocated(rectangles)) then
            do j_element = size(rectangles, 2), 1, -1
                do i_element = size(rectangles, 1), 1, -1
                    call rectangles(i_element, j_element)%writer%destroy()
                end do
            end do
            deallocate(rectangles)
        end if
    end subroutine destroy_rectangle

    subroutine create_transmutations_successes(successes, paths, filename, components)
        class(Abstract_Rectangle_Writer), allocatable, intent(out) :: successes(:)
        type(String_Wrapper), intent(in) :: paths(:)
        character(len=*), intent(in) :: filename
        type(Component_Wrapper), intent(in) :: components(:, :)

        logical :: selectors(size(components, 1), size(components, 1), size(components, 2))
        logical :: some_couples_can_exchange(size(selectors, 3)), can_exchange_ij
        integer :: i_box, i_component, j_component

        some_couples_can_exchange = .false.
        do i_box = 1, size(selectors, 3)
            do j_component = 1, size(selectors, 2)
            do i_component = 1, size(selectors, 1)
                can_exchange_ij = &
                    component_can_exchange(components(i_component, i_box)%chemical_potential) .and.&
                    component_can_exchange(components(j_component, i_box)%chemical_potential) .and.&
                    i_component /= j_component
                some_couples_can_exchange(i_box) = some_couples_can_exchange(i_box) .or. &
                    can_exchange_ij
                selectors(i_component, j_component, i_box) = can_exchange_ij
            end do
            end do
        end do

        if (all(some_couples_can_exchange)) then
            allocate(Concrete_Rectangle_Writer :: successes(size(paths)))
        else
            allocate(Null_Rectangle_Writer :: successes(size(paths)))
        end if

        do i_box = 1, size(successes)
            call successes(i_box)%construct(paths(i_box)%string//filename, selectors(:, :, i_box))
        end do
    end subroutine create_transmutations_successes

    subroutine destroy_line(rectangles)
        class(Abstract_Rectangle_Writer), allocatable, intent(out) :: rectangles(:)

        integer :: i_element

        if (allocated(rectangles)) then
            do i_element = size(rectangles), 1, -1
                call rectangles(i_element)%destroy()
            end do
            deallocate(rectangles)
        end if
    end subroutine destroy_line

end module procedures_rectangle_writer_factory
