module types_changes_wrapper

use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_random_coordinates, only: Abstract_Random_Coordinates
use classes_coordinates_copier, only: Abstract_Coordinates_Copier
use types_changes_component_wrapper, only: Changes_Component_Wrapper

implicit none

private

    type, public :: Changes_Wrapper
        class(Abstract_Parallelepiped_Domain), allocatable :: exchange_domain
        class(Abstract_Random_Coordinates), allocatable :: random_position, random_orientation
        class(Abstract_Coordinates_Copier), allocatable :: position_copier, orientation_copier
        type(Changes_Component_Wrapper), allocatable :: components(:)
    end type Changes_Wrapper

end module types_changes_wrapper
