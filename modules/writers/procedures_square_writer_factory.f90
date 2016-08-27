module procedures_square_writer_factory

use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only: component_can_exchange
use classes_square_writer, only: Abstract_Square_Writer, Concrete_Square_Writer, Null_Square_Writer

implicit none

private
public :: create_transmutations, destroy

contains

    subroutine create_transmutations(transmutations, components, filename)
        class(Abstract_Square_Writer), allocatable, intent(out) :: transmutations
        type(Component_Wrapper), intent(in) :: components(:)
        character(len=*), intent(in) :: filename

        logical :: selectors(size(components), size(components))
        logical :: some_couples_can_exchange, can_exchange_ij
        integer :: i_component, j_component

        some_couples_can_exchange = .false.
        do j_component = 1, size(selectors, 2)
            do i_component = 1, size(selectors, 1)
                can_exchange_ij = &
                    component_can_exchange(components(i_component)%chemical_potential) .and. &
                    component_can_exchange(components(j_component)%chemical_potential) .and. &
                    i_component /= j_component
                some_couples_can_exchange = some_couples_can_exchange .or. can_exchange_ij
                selectors(i_component, j_component) = can_exchange_ij
            end do
        end do
        if (some_couples_can_exchange) then
            allocate(Concrete_Square_Writer :: transmutations)
        else
            allocate(Null_Square_Writer :: transmutations)
        end if
        call transmutations%construct(selectors, filename)
    end subroutine create_transmutations

    subroutine destroy(square)
        class(Abstract_Square_Writer), allocatable, intent(out) :: square

        if (allocated(square)) then
            call square%destroy()
            deallocate(square)
        end if
    end subroutine destroy

end module procedures_square_writer_factory
