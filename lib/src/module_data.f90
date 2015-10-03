module module_data

use procedures_errors, only: error_exit

implicit none

private
public :: test_file_exists, test_data_found, test_empty_string

contains

    subroutine test_file_exists(filename)
        character(len=*), intent(in) :: filename

        logical :: file_exists

        inquire(file=filename, exist=file_exists)
        if (.not.file_exists) then
            call error_exit(filename//" doesn't exist.")
        end if
    end subroutine test_file_exists

    subroutine test_data_found(field, found)
        character(len=*), intent(in) :: field
        logical, intent(in) :: found

        if (.not.found) then
            call error_exit(trim(field)//" not found.")
        end if
    end subroutine test_data_found

    subroutine test_empty_string(field, string)
        character(len=*), intent(in) :: field
        character(len=*) :: string

        if (len(string) == 0) then
            call error_exit(trim(field)//": string is empty.")
        end if
    end subroutine test_empty_string

end module module_data
