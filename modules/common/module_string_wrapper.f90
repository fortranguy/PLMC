module module_string_wrapper

use classes_number_to_string, only: Abstract_Number_to_String

implicit none

private
public :: strings_wrapper_destroy

    type, public :: String_Wrapper
        class(Abstract_Number_to_String), allocatable :: string
    end type String_Wrapper

    type, public :: Strings_Line
        type(String_Wrapper), allocatable :: line(:)
    end type Strings_Line

contains

    subroutine strings_wrapper_destroy(strings)
        type(String_Wrapper), allocatable, intent(inout) :: strings(:)

        integer :: i_string

        if (allocated(strings)) then
            do i_string = size(strings), 1, -1
                if (allocated(strings(i_string)%string)) deallocate(strings(i_string)%string)
            end do
        end if
    end subroutine strings_wrapper_destroy

end module module_string_wrapper
