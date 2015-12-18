module class_number_to_string

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: max_word_length

implicit none

private

    type, abstract, public :: Abstract_Number_to_String
    contains
        generic :: get => get_real_scalar, get_real_array, get_integer
        procedure, nopass, private :: get_real_scalar => Abstract_get_real_scalar
        procedure, nopass, private :: get_real_array => Abstract_get_real_array
        procedure, nopass, private :: get_integer => Abstract_get_integer
    end type Abstract_Number_to_String

    type, extends(Abstract_Number_to_String), public :: Concrete_Number_to_String

    end type Concrete_Number_to_String

    type, extends(Abstract_Number_to_String), public :: Null_Number_to_String
    contains
        procedure, nopass, private :: get_real_scalar => Null_get_real_scalar
        procedure, nopass, private :: get_real_array => Null_get_real_array
        procedure, nopass, private :: get_integer => Null_get_integer
    end type Null_Number_to_String

contains

!implementation Abstract_Number_to_String

    function Abstract_get_real_scalar(number) result(string)
        character(len=:), allocatable :: string
        real(DP), intent(in) :: number

        character(len=max_word_length) :: big_string

        write(big_string, *) number
        string = trim(big_string)
    end function Abstract_get_real_scalar
    !elemental: works with gfortran 4.9 / doesn't work with ifort 13

    function Abstract_get_real_array(number) result(string)
        character(len=:), allocatable :: string
        real(DP), intent(in) :: number(:)

        character(len=max_word_length) :: big_string

        write(big_string, *) number
        string = trim(big_string)
    end function Abstract_get_real_array

    function Abstract_get_integer(number) result(string)
        character(len=:), allocatable :: string
        integer, intent(in) :: number

        character(len=max_word_length) :: big_string

        write(big_string, *) number
        string = trim(adjustl(big_string))
    end function Abstract_get_integer

!end implementation Abstract_Number_to_String

!implementation Null_Number_to_String

    function Null_get_real_scalar(number) result(string)
        character(len=:), allocatable :: string
        real(DP), intent(in) :: number
        string = ""
    end function Null_get_real_scalar

    function Null_get_real_array(number) result(string)
        character(len=:), allocatable :: string
        real(DP), intent(in) :: number(:)
        string = ""
    end function Null_get_real_array

    function Null_get_integer(number) result(string)
        character(len=:), allocatable :: string
        integer, intent(in) :: number
        string = ""
    end function Null_get_integer

!end implementation Null_Number_to_String

end module class_number_to_string
