module procedures_checks

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use procedures_errors, only: warning_continue, error_exit

implicit none

private
public check_positive_real

interface check_positive_real
    module procedure check_positive_real_scalar
    module procedure check_positive_real_array
end interface check_positive_real

contains

    subroutine check_positive_real_scalar(class_name, real_name, real_scalar)
        character(len=*), intent(in) :: class_name, real_name
        real(DP), intent(in) :: real_scalar

        if (real_scalar < 0._DP) then
            call error_exit(class_name//": "//real_name//" is negative.")
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
    
end module procedures_checks
