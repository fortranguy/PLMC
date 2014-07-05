module module_data

use, intrinsic :: iso_fortran_env, only: error_unit
use json_module, only: json_file

implicit none

private
public json_get_string, test_data_found

contains

    subroutine json_get_string(json, data_name, string)
    
        type(json_file), intent(inout) :: json
        character(len=*), intent(in) :: data_name
        character(len=:), allocatable, intent(out) :: string
        
        character(len=:), allocatable :: string_dummy
        logical :: found
        
        call json%get(data_name, string_dummy, found)
        call test_data_found(data_name, found)
        
        write(*, *) "string_dummy: ", string_dummy
        
        if (len(string_dummy) == 0) then
            write(error_unit, *) "String is empty."
            error stop
        else
            string = string_dummy
            if (allocated(string_dummy)) deallocate(string_dummy)
        end if
        
        string = "blog !"
        
    end subroutine json_get_string 

    subroutine test_data_found(data_name, found)
    
        character(len=*), intent(in) :: data_name
        logical, intent(in) :: found

        if (.not.found) then
            write(error_unit, *) trim(data_name), " not found."
            error stop
        end if
        
    end subroutine test_data_found

end module module_data
