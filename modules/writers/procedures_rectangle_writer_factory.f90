module procedures_rectangle_writer_factory

use types_string_wrapper, only: String_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only: component_can_exchange
use classes_rectangle_writer, only: Abstract_Rectangle_Writer, Concrete_Rectangle_Writer, &
    Null_Rectangle_Writer

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_transmutations_successes
end interface create

interface destroy
    module procedure :: destroy_line
end interface destroy

contains

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
