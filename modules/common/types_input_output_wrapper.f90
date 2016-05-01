module types_input_output_wrapper

use json_module, only: json_file
use types_readers_wrapper, only: Readers_Wrapper
use types_writers_wrapper, only: Writers_Wrapper

private

    type, public :: Input_Output_Wrapper
        type(json_file) :: input_data
        type(Readers_Wrapper) :: readers
        type(Writers_Wrapper) :: writers
    end type Input_Output_Wrapper

end module types_input_output_wrapper
