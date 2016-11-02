module procedures_hetero_couples_factory

use classes_hetero_couples, only: Abstract_Hetero_Couples, Half_Hetero_Couples, &
    Full_Hetero_Couples, Null_Hetero_Couples

implicit none

private
public :: create_half, create_full, destroy

interface create_full
    module procedure :: create_full_line
    module procedure :: create_full_element
end interface create_full

interface create_half
    module procedure :: create_half_line
    module procedure :: create_half_element
end interface create_half

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_line
end interface destroy

contains

    subroutine create_half_line(couples, num_elements, num_partners)
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
    end subroutine create_half_line

    subroutine create_half_element(couples, num_partners)
        class(Abstract_Hetero_Couples), allocatable, intent(out) :: couples
        integer, intent(in) :: num_partners

        if (num_partners > 1) then
            allocate(Half_Hetero_Couples :: couples)
        else
            allocate(Null_Hetero_Couples :: couples)
        end if

        call couples%construct(num_partners)
    end subroutine create_half_element

    subroutine create_full_line(couples, num_elements, num_partners)
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
    end subroutine create_full_line

    subroutine create_full_element(couples, num_partners)
        class(Abstract_Hetero_Couples), allocatable, intent(out) :: couples
        integer, intent(in) :: num_partners

        if (num_partners > 1) then
            allocate(Full_Hetero_Couples :: couples)
        else
            allocate(Null_Hetero_Couples :: couples)
        end if

        call couples%construct(num_partners)
    end subroutine create_full_element

    subroutine destroy_line(couples)
        class(Abstract_Hetero_Couples), allocatable, intent(inout) :: couples(:)

        integer :: i_element

        if (allocated(couples)) then
            do i_element = size(couples), 1, -1
                call couples(i_element)%destroy()
            end do
            deallocate(couples)
        end if
    end subroutine destroy_line

    subroutine destroy_element(couples)
        class(Abstract_Hetero_Couples), allocatable, intent(inout) :: couples

        if (allocated(couples)) then
            call couples%destroy()
            deallocate(couples)
        end if
    end subroutine destroy_element

end module procedures_hetero_couples_factory
