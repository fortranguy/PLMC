module procedures_string_factory

use types_string_wrapper, only: String_Wrapper, Number_to_String_Wrapper, Number_to_String_Line

implicit none

private
public :: destroy

interface destroy
    module procedure :: destroy_number_triangle
    module procedure :: destroy_number_line
    module procedure :: destroy_line
end interface destroy

contains

    subroutine destroy_number_triangle(strings)
        type(Number_to_String_Line), allocatable, intent(inout) :: strings(:)

        integer :: i_string

        if (allocated(strings)) then
            do i_string = size(strings), 1, -1
                call destroy(strings(i_string)%line)
            end do
            deallocate(strings)
        end if
    end subroutine destroy_number_triangle

    subroutine destroy_number_line(strings)
        type(Number_to_String_Wrapper), allocatable, intent(inout) :: strings(:)

        integer :: i_string

        if (allocated(strings)) then
            do i_string = size(strings), 1, -1
                if (allocated(strings(i_string)%string)) deallocate(strings(i_string)%string)
            end do
        end if
    end subroutine destroy_number_line

    !> @note Fortran needs generics.
    subroutine destroy_line(strings)
        type(String_Wrapper), allocatable, intent(inout) :: strings(:)

        integer :: i_string

        if (allocated(strings)) then
            do i_string = size(strings), 1, -1
                if (allocated(strings(i_string)%string)) deallocate(strings(i_string)%string)
            end do
        end if
    end subroutine destroy_line

end module procedures_string_factory
