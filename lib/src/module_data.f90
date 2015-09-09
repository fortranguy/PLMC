module module_data

use procedures_errors, only: error_exit

implicit none

private
public data_filename, data_post_filename, report_filename, report_post_filename, &
       test_file_exists, test_data_found, test_empty_string

    character(len=*), parameter :: data_filename = "data.json"
    character(len=*), parameter :: data_post_filename = "data_post.json"
    character(len=*), parameter :: report_filename = "report.json"
    character(len=*), parameter :: report_post_filename = "report_post.json"

contains

    subroutine test_file_exists(filename)
        character(len=*), intent(in) :: filename
        
        logical :: file_exists
        
        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            call error_exit(filename//" doesn't exist.")
        end if
    end subroutine test_file_exists

    subroutine test_data_found(data_field, found)
        character(len=*), intent(in) :: data_field
        logical, intent(in) :: found

        if (.not. found) then
            call error_exit(trim(data_field)//" not found.")
        end if
    end subroutine test_data_found
    
    subroutine test_empty_string(data_field, string)
        character(len=*), intent(in) :: data_field
        character(len=*) :: string
        
        if (len(string) == 0) then
            call error_exit(trim(data_field)//": string is empty.")
        end if
    end subroutine test_empty_string

end module module_data
