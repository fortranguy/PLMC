module module_string_wrapper

use class_number_to_string, only: Abstract_Number_to_String

implicit none

private
public :: strings_wrapper_destroy

    type, public :: String_Wrapper
        class(Abstract_Number_to_String), allocatable :: string
    end type String_Wrapper

contains

    subroutine strings_wrapper_destroy(strings)
        type(String_Wrapper), allocatable, intent(inout) :: strings(:)

        integer :: i_string

        if (allocated(strings)) then
            do i_string = 1, size(strings)
                if (allocated(strings(i_string)%string)) then
                    deallocate(strings(i_string)%string)
                end if
            end do
        end if
    end subroutine strings_wrapper_destroy

end module module_string_wrapper
