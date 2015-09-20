module procedures_checks

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use data_geometry, only: num_dimensions
use procedures_errors, only: warning_continue, error_exit

implicit none

private
public :: check_in_range, check_3d_array, check_positive, check_norm, check_adaptation_factor

interface check_3d_array
    module procedure check_integer_3d_array
    module procedure check_real_3d_array
end interface check_3d_array

interface check_positive
    module procedure check_positive_integer_scalar
    module procedure check_positive_integer_array
    module procedure check_positive_real_scalar
    module procedure check_positive_real_array
end interface check_positive

contains

    subroutine check_in_range(class_name, integer_max, integer_name, integer_value)
        character(len=*), intent(in) :: class_name, integer_name
        integer, intent(in) :: integer_max, integer_value

        character(len=1024) :: string_i

        write(string_i, *) integer_value
        if (integer_value < 1 .or. integer_max < integer_value) then
            call error_exit(class_name//": "//integer_name//"="//trim(adjustl(string_i))//&
                            " is out of range.")
        end if
    end subroutine check_in_range

!implementation check_3d_array

    subroutine check_integer_3d_array(class_name, integer_name, integer_array)
        character(len=*), intent(in) :: class_name, integer_name
        integer, intent(in) :: integer_array(:)

        if (size(integer_array) /= num_dimensions) then
            call error_exit(class_name//": "//integer_name//" has wrong number of dimensions (size).")
        end if
    end subroutine check_integer_3d_array

    subroutine check_real_3d_array(class_name, real_name, real_array)
        character(len=*), intent(in) :: class_name, real_name
        real(DP), intent(in) :: real_array(:)

        if (size(real_array) /= num_dimensions) then
            call error_exit(class_name//": "//real_name//" has wrong number of dimensions (size).")
        end if
    end subroutine check_real_3d_array

!end implementation check_3d_array

!implementation check_positive

    subroutine check_positive_integer_scalar(class_name, integer_name, integer_scalar)
        character(len=*), intent(in) :: class_name, integer_name
        integer, intent(in) :: integer_scalar

        character(len=1024) :: string_i

        write(string_i, *) integer_scalar
        if (integer_scalar < 0) then
            call error_exit(class_name//": "//integer_name//"="//trim(adjustl(string_i))//&
                            " is negative.")
        end if
        if (integer_scalar == 0) then
            call warning_continue(class_name//": "//integer_name//" is zero.")
        end if
    end subroutine check_positive_integer_scalar

    subroutine check_positive_integer_array(class_name, integer_name, integer_array)
        character(len=*), intent(in) :: class_name, integer_name
        integer, intent(in) :: integer_array(:)

        integer :: i_dimension
        character(len=1024) :: string_i

        do i_dimension = 1, size(integer_array)
            write(string_i, *) i_dimension
            call check_positive_integer_scalar(class_name, &
                                               integer_name//"("//trim(adjustl(string_i))//")", &
                                               integer_array(i_dimension))
        end do
    end subroutine check_positive_integer_array

    subroutine check_positive_real_scalar(class_name, real_name, real_scalar)
        character(len=*), intent(in) :: class_name, real_name
        real(DP), intent(in) :: real_scalar

        character(len=1024) :: string_real

        write(string_real, *) real_scalar

        if (real_scalar < 0._DP) then
            call error_exit(class_name//": "//real_name//"="//trim(adjustl(string_real))//&
                            " is negative.")
        end if
        if (real_scalar < real_zero) then
            call warning_continue(class_name//": "//real_name//" may be too small.")
        end if
    end subroutine check_positive_real_scalar

    subroutine check_positive_real_array(class_name, real_name, real_array)
        character(len=*), intent(in) :: class_name, real_name
        real(DP), intent(in) :: real_array(:)

        integer :: i_dimension
        character(len=1024) :: string_i

        do i_dimension = 1, size(real_array)
            write(string_i, *) i_dimension
            call check_positive_real_scalar(class_name, &
                                            real_name//"("//trim(adjustl(string_i))//")", &
                                            real_array(i_dimension))
        end do
    end subroutine check_positive_real_array

!end implementation check_positive

    subroutine check_norm(class_name, vector_name, vector)
        character(len=*), intent(in) :: class_name, vector_name
        real(DP), intent(in) :: vector(:)

        character(len=1024) :: string_vector

        if (abs(norm2(vector) - 1.0_DP) > real_zero) then
            write(string_vector, *) norm2(vector)
            call warning_continue(class_name//": "//vector_name//" may not be normed "//&
            "("//trim(adjustl(string_vector))//").")
        end if
    end subroutine check_norm

    subroutine check_adaptation_factor(class_name, value)
        character(len=*), intent(in) :: class_name
        real(DP), intent(in) :: value

        if (value < 1.0_DP) then
            call error_exit(class_name//": adaptation_factor is less than 1.0.")
        end if
        if (abs(value - 1.0_DP) < real_zero) then
            call warning_continue(class_name//": adaptation_factor is 1.0.")
        end if
    end subroutine check_adaptation_factor

end module procedures_checks
