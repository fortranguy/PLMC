module module_data

use, intrinsic :: iso_fortran_env, only: error_unit

implicit none

private
public test_data_file_exists, test_data_found, test_empty_string

contains

    subroutine test_data_file_exists(data_filename)

        character(len=*), intent(in) :: data_filename
        
        logical :: data_exists
        
        inquire(file=data_filename, exist=data_exists)
        
        if (.not. data_exists) then
            write(error_unit, *) data_filename, " doesn't exist."
            error stop
        end if

    end subroutine test_data_file_exists

    subroutine test_data_found(data_name, found)
    
        character(len=*), intent(in) :: data_name
        logical, intent(in) :: found

        if (.not.found) then
            write(error_unit, *) trim(data_name), " not found."
            error stop
        end if
        
    end subroutine test_data_found
    
    subroutine test_empty_string(data_name, string)
    
        character(len=*), intent(in) :: data_name
        character(len=*) :: string
        
        if (len(string) == 0) then
            write(error_unit, *) trim(data_name), ": string is empty."
            error stop
        end if
        
    end subroutine test_empty_string

end module module_data
