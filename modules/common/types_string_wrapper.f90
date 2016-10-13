module types_string_wrapper

use classes_number_to_string, only: Abstract_Number_to_String

implicit none

private

    type, public :: String_Wrapper
        character(len=:), allocatable :: string
    end type String_Wrapper

    type, public :: Number_to_String_Wrapper
        class(Abstract_Number_to_String), allocatable :: string
    end type Number_to_String_Wrapper

    type, public :: Number_to_String_Line
        type(Number_to_String_Wrapper), allocatable :: line(:)
    end type Number_to_String_Line

end module types_string_wrapper
