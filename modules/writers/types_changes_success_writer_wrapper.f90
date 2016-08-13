module types_changes_success_writer_wrapper

use classes_changes_success_writer, only: Abstract_Changes_Success_Writer

implicit none

private

    type, public :: Changes_Success_Writer_Wrapper
        class(Abstract_Changes_Success_Writer), allocatable :: writer
    end type Changes_Success_Writer_Wrapper

end module types_changes_success_writer_wrapper
