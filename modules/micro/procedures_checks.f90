module procedures_checks

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, max_word_length, real_zero
use procedures_errors, only: warning_continue, error_exit
use types_potential_domain, only: Concrete_Potential_Domain

implicit none

private
public :: check_file_exists, check_data_found, check_string_not_empty, &
    check_in_range, check_3d_array, check_positive, check_norm, check_increase_factor, &
    check_potential_domain, check_ratio

interface check_3d_array
    module procedure :: check_integer_3d_array
    module procedure :: check_real_3d_array
end interface check_3d_array

interface check_positive
    module procedure :: check_positive_integer_scalar
    module procedure :: check_positive_integer_array
    module procedure :: check_positive_real_scalar
    module procedure :: check_positive_real_array
end interface check_positive

contains

    subroutine check_file_exists(filename)
        character(len=*), intent(in) :: filename

        logical :: file_exists

        inquire(file=filename, exist=file_exists)
        if (.not.file_exists) then
            call error_exit(filename//" doesn't exist.")
        end if
    end subroutine check_file_exists

    subroutine check_data_found(field, found)
        character(len=*), intent(in) :: field
        logical, intent(in) :: found

        if (.not.found) then
            call error_exit(trim(field)//" not found.")
        end if
    end subroutine check_data_found

    subroutine check_string_not_empty(field, string)
        character(len=*), intent(in) :: field
        character(len=*) :: string

        if (len(string) == 0) then
            call error_exit(trim(field)//": string is empty.")
        end if
    end subroutine check_string_not_empty

    subroutine check_in_range(context, integer_max, integer_name, integer_value)
        character(len=*), intent(in) :: context, integer_name
        integer, intent(in) :: integer_max, integer_value

        character(len=max_word_length) :: string_i

        write(string_i, *) integer_value
        if (integer_value < 1 .or. integer_max < integer_value) then
            call error_exit(context//": "//integer_name//"="//trim(adjustl(string_i))//&
                            " is out of range.")
        end if
    end subroutine check_in_range

!implementation check_3d_array

    subroutine check_integer_3d_array(context, integer_name, integer_array)
        character(len=*), intent(in) :: context, integer_name
        integer, intent(in) :: integer_array(:)

        if (size(integer_array) /= num_dimensions) then
            call error_exit(context//": "//integer_name//" has wrong number of dimensions (size).")
        end if
    end subroutine check_integer_3d_array

    subroutine check_real_3d_array(context, real_name, real_array)
        character(len=*), intent(in) :: context, real_name
        real(DP), intent(in) :: real_array(:)

        if (size(real_array) /= num_dimensions) then
            call error_exit(context//": "//real_name//" has wrong number of dimensions (size).")
        end if
    end subroutine check_real_3d_array

!end implementation check_3d_array

!implementation check_positive

    subroutine check_positive_integer_scalar(context, integer_name, integer_scalar)
        character(len=*), intent(in) :: context, integer_name
        integer, intent(in) :: integer_scalar

        character(len=max_word_length) :: string_i

        write(string_i, *) integer_scalar
        if (integer_scalar < 0) then
            call error_exit(context//": "//integer_name//"="//trim(adjustl(string_i))//&
                            " is negative.")
        end if
        if (integer_scalar == 0) then
            call warning_continue(context//": "//integer_name//" is zero.")
        end if
    end subroutine check_positive_integer_scalar

    subroutine check_positive_integer_array(context, integer_name, integer_array)
        character(len=*), intent(in) :: context, integer_name
        integer, intent(in) :: integer_array(:)

        integer :: i_dimension
        character(len=max_word_length) :: string_i

        do i_dimension = 1, size(integer_array)
            write(string_i, *) i_dimension
            call check_positive_integer_scalar(context, &
                                               integer_name//"("//trim(adjustl(string_i))//")", &
                                               integer_array(i_dimension))
        end do
    end subroutine check_positive_integer_array

    subroutine check_positive_real_scalar(context, real_name, real_scalar)
        character(len=*), intent(in) :: context, real_name
        real(DP), intent(in) :: real_scalar

        character(len=max_word_length) :: string_real

        write(string_real, *) real_scalar

        if (real_scalar < 0._DP) then
            call error_exit(context//": "//real_name//"="//trim(adjustl(string_real))//&
                            " is negative.")
        end if
        if (real_scalar < real_zero) then
            call warning_continue(context//": "//real_name//" may be too small.")
        end if
    end subroutine check_positive_real_scalar

    subroutine check_positive_real_array(context, real_name, real_array)
        character(len=*), intent(in) :: context, real_name
        real(DP), intent(in) :: real_array(:)

        integer :: i_dimension
        character(len=max_word_length) :: string_i

        do i_dimension = 1, size(real_array)
            write(string_i, *) i_dimension
            call check_positive_real_scalar(context, &
                                            real_name//"("//trim(adjustl(string_i))//")", &
                                            real_array(i_dimension))
        end do
    end subroutine check_positive_real_array

!end implementation check_positive

    subroutine check_norm(context, vector_name, vector)
        character(len=*), intent(in) :: context, vector_name
        real(DP), intent(in) :: vector(:)

        character(len=max_word_length) :: string_vector

        if (abs(norm2(vector) - 1.0_DP) > real_zero) then
            write(string_vector, *) norm2(vector)
            call warning_continue(context//": "//vector_name//" may not be normed "//&
            "("//trim(adjustl(string_vector))//").")
        end if
    end subroutine check_norm

    subroutine check_increase_factor(context, increase_factor_name, increase_factor)
        character(len=*), intent(in) :: context, increase_factor_name
        real(DP), intent(in) :: increase_factor

        call check_positive(context, increase_factor_name, increase_factor)
        if (increase_factor < 1.0_DP) then
            call error_exit(context//": "//increase_factor_name//" is less than 1.0.")
        end if
        if (abs(increase_factor - 1.0_DP) < real_zero) then
            call warning_continue(context//": "//increase_factor_name//" is 1.0.")
        end if
    end subroutine check_increase_factor

    subroutine check_potential_domain(context, domain_name, domain)
        character(len=*), intent(in) :: context, domain_name
        type(Concrete_Potential_Domain), intent(in) :: domain

        real(DP) :: distance_range
        call check_positive(context, domain_name//"%min", domain%min)
        call check_positive(context, domain_name//"%max", domain%max)
        if (domain%min > domain%max) then
            call error_exit(context//": "//domain_name//"%min > domain%max.")
        end if
        call check_positive(context, "domain%delta", domain%delta)
        distance_range = domain%max - domain%min
        if (distance_range < real_zero) then
            call warning_continue(context//"distance_range may be too small.")
        end if
        if (distance_range / domain%delta < 1._DP) then
            call warning_continue(context//domain_name//"%delta may be too big.")
        end if
    end subroutine check_potential_domain

    subroutine check_ratio(context, ratio_name, ratio)
        character(len=*), intent(in) :: context, ratio_name
        real(DP), intent(in) :: ratio

        if (ratio < 0._DP .or. 1._DP < ratio) then
            call error_exit(context//": "//ratio_name//" must be between 0.0 and 1.0.")
        end if
    end subroutine check_ratio

end module procedures_checks