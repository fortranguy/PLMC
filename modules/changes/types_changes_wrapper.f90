module types_changes_wrapper

use classes_random_coordinates, only: Abstract_Random_Coordinates
use types_changes_component_wrapper, only: Changes_Component_Wrapper

implicit none

    type, public :: Changes_Wrapper
        class(Abstract_Random_Coordinates), allocatable :: random_position, random_orientation
        type(Changes_Component_Wrapper), allocatable :: components(:)
    end type Changes_Wrapper

end module types_changes_wrapper
