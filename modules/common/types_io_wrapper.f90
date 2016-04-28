module types_io_wrapper

use types_readers_wrapper, only: Readers_Wrapper
use types_writers_wrapper, only: Writers_Wrapper

private

    type, public :: IO_Wrapper
        type(Readers_Wrapper) :: readers
        type(Writers_Wrapper) :: writers
    end type IO_Wrapper

end module types_io_wrapper
