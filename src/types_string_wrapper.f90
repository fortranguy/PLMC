module types_string_wrapper

use classes_number_to_string, only: Abstract_Number_to_String

implicit none

private

    type, public :: String_Wrapper
        character(len=:), allocatable :: string
    end type String_Wrapper

end module types_string_wrapper
