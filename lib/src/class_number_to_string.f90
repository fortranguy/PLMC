module class_number_to_string

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, abstract, public :: Abstract_Number_to_String
    contains
        generic :: get => get_real, get_integer
        procedure, nopass, private :: get_real => Abstract_Number_to_String_get_real
        procedure, nopass, private :: get_integer => Abstract_Number_to_String_get_integer
    end type Abstract_Number_to_String

    type, extends(Abstract_Number_to_String), public :: Concrete_Number_to_String

    end type Concrete_Number_to_String

    type, extends(Abstract_Number_to_String), public :: Null_Number_to_String
    contains
        procedure, nopass, private :: get_real => Null_Number_to_String_get_real
        procedure, nopass, private :: get_integer => Null_Number_to_String_get_integer
    end type Null_Number_to_String

contains

!implementation Abstract_Number_to_String

    function Abstract_Number_to_String_get_real(number) result(string)
        character(len=:), allocatable :: string
        real(DP), intent(in) :: number

        character(len=1024) :: big_string

        write(big_string, *) number
        string = trim(big_string)
    end function Abstract_Number_to_String_get_real

    function Abstract_Number_to_String_get_integer(number) result(string)
        character(len=:), allocatable :: string
        integer, intent(in) :: number

        character(len=1024) :: big_string

        write(big_string, *) number
        string = trim(big_string)
    end function Abstract_Number_to_String_get_integer

!end implementation Abstract_Number_to_String

!implementation Null_Number_to_String

    function Null_Number_to_String_get_real(number) result(string)
        character(len=:), allocatable :: string
        real(DP), intent(in) :: number
        string = ""
    end function Null_Number_to_String_get_real

    function Null_Number_to_String_get_integer(number) result(string)
        character(len=:), allocatable :: string
        integer, intent(in) :: number
        string = ""
    end function Null_Number_to_String_get_integer

!end implementation Null_Number_to_String

end module class_number_to_string
