module module_data

use, intrinsic :: iso_fortran_env, only: error_unit

implicit none

private
public data_filename, data_post_filename, report_filename, report_post_filename, &
       spheres1_object_field, spheres2_object_field, &
       test_file_exists, test_data_found, test_empty_string

    character(len=*), parameter :: data_filename = "data.json"
    character(len=*), parameter :: data_post_filename = "data_post.json"
    character(len=*), parameter :: report_filename = "report.json"
    character(len=*), parameter :: report_post_filename = "report_post.json"
    
    character(len=*), parameter :: spheres1_object_field = "Spheres 1"
    character(len=*), parameter :: spheres2_object_field = "Spheres 2"

contains

    subroutine test_file_exists(filename)
        character(len=*), intent(in) :: filename
        
        logical :: data_exists
        
        inquire(file=filename, exist=data_exists)
        
        if (.not. data_exists) then
            write(error_unit, *) filename, " doesn't exist."
            error stop
        end if
    end subroutine test_file_exists

    subroutine test_data_found(data_field, found)
        character(len=*), intent(in) :: data_field
        logical, intent(in) :: found

        if (.not.found) then
            write(error_unit, *) trim(data_field), " not found."
            error stop
        end if        
    end subroutine test_data_found
    
    subroutine test_empty_string(data_field, string)
        character(len=*), intent(in) :: data_field
        character(len=*) :: string
        
        if (len(string) == 0) then
            write(error_unit, *) trim(data_field), ": string is empty."
            error stop
        end if        
    end subroutine test_empty_string

end module module_data
