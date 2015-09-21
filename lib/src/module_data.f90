module module_data

use procedures_errors, only: error_exit
use json_module, only: json_file

implicit none

private
public :: data_filename, data_post_filename, report_filename, report_post_filename, &
    test_file_exists, test_data_found, test_empty_string

    character(len=*), parameter :: data_filename = "data.json"
    character(len=*), parameter :: data_post_filename = "data_post.json"
    character(len=*), parameter :: report_filename = "report.json"
    character(len=*), parameter :: report_post_filename = "report_post.json"

    type, public :: Concrete_Input_Data
        type(json_file) :: json
        character(len=:), allocatable :: prefix
    end type Concrete_Input_Data

contains

    subroutine test_file_exists(filename)
        character(len=*), intent(in) :: filename

        logical :: file_exists

        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            call error_exit(filename//" doesn't exist.")
        end if
    end subroutine test_file_exists

    subroutine test_data_found(field, found)
        character(len=*), intent(in) :: field
        logical, intent(in) :: found

        if (.not. found) then
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
