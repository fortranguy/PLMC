module procedures_hetero_couples_factory

use classes_hetero_couples, only: Abstract_Hetero_Couples, Half_Hetero_Couples, &
    Full_Hetero_Couples, Null_Hetero_Couples

implicit none

private
public :: create_half, create_full, destroy

contains

    subroutine create_half(couples, num_elements, num_partners)
        class(Abstract_Hetero_Couples), allocatable, intent(out) :: couples(:)
        integer, intent(in) :: num_elements, num_partners

        integer :: i_element

        if (num_partners > 1) then
            allocate(Half_Hetero_Couples :: couples(num_elements))
        else
            allocate(Null_Hetero_Couples :: couples(num_elements))
        end if

        do i_element = 1, size(couples)
            call couples(i_element)%construct(num_partners)
        end do
    end subroutine create_half

    subroutine create_full(couples, num_elements, num_partners)
        class(Abstract_Hetero_Couples), allocatable, intent(out) :: couples(:)
        integer, intent(in) :: num_elements, num_partners

        integer :: i_element

        if (num_partners > 1) then
            allocate(Full_Hetero_Couples :: couples(num_elements))
        else
            allocate(Null_Hetero_Couples :: couples(num_elements))
        end if

        do i_element = 1, size(couples)
            call couples(i_element)%construct(num_partners)
        end do
    end subroutine create_full

    subroutine destroy(couples)
        class(Abstract_Hetero_Couples), allocatable, intent(inout) :: couples(:)

        integer :: i_element

        if (allocated(couples)) then
            do i_element = size(couples), 1, -1
                call couples(i_element)%destroy()
            end do
            deallocate(couples)
        end if
    end subroutine destroy

end module procedures_hetero_couples_factory
