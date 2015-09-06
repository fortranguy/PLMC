module module_error

use, intrinsic :: iso_fortran_env, only: error_unit

private
public warning_continue, error_exit

contains

    subroutine warning_continue(message)
        character(len=*), intent(in) :: message

        write(error_unit, *) "Warning: ", message
    end subroutine warning_continue

    subroutine error_exit(message)
        character(len=*), intent(in) :: message
        
        write(error_unit, *) "Fatal Error: ", message
        error stop
    end subroutine error_exit

end module module_error