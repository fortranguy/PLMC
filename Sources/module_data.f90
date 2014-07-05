module module_data

use, intrinsic :: iso_fortran_env, only: error_unit
use json_module, only: json_file

implicit none

private
public json_get_string, test_data_found

contains

    function json_get_string(json, data_name)
    
        type(json_file), intent(inout) :: json
        character(len=*), intent(in) :: data_name
        character(len=:), allocatable :: json_get_string
        
        logical :: found
        
        call json%get(data_name, json_get_string, found)
        call test_data_found(data_name, found)
        
        if (len(json_get_string) == 0) then
            write(error_unit, *) trim(data_name), ": String is empty."
            error stop
        end if
        
    end function json_get_string 

    subroutine test_data_found(data_name, found)
    
        character(len=*), intent(in) :: data_name
        logical, intent(in) :: found

        if (.not.found) then
            write(error_unit, *) trim(data_name), " not found."
            error stop
        end if
        
    end subroutine test_data_found

end module module_data
