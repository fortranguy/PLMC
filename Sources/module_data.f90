module module_data

use, intrinsic :: iso_fortran_env, only: error_unit

implicit none

private
public test_data_not_found

contains

    subroutine test_data_not_found(data_name, found)
        character(len=*), intent(in) :: data_name
        logical, intent(in) :: found

        if (.not.found) then
            write(error_unit, *) trim(data_name), " not found."
            error stop
        end if
        
    end subroutine test_data_not_found

end module module_data
