module types_changes_wrapper

use types_changes_component_wrapper, only: Changes_Component_Wrapper

implicit none

    type, public :: Changes_Wrapper
        type(Changes_Component_Wrapper), allocatable :: components(:)
    end type Changes_Wrapper

end module types_changes_wrapper
