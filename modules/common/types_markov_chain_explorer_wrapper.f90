module types_markov_chain_explorer_wrapper

use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_random_coordinates, only: Abstract_Random_Coordinates
use classes_widom_method, only: Abstract_Widom_Method

implicit none

private

    type, public :: Markov_Chain_Explorer_Wrapper
        class(Abstract_Parallelepiped_Domain), allocatable :: widom_domain
        class(Abstract_Random_Coordinates), allocatable :: random_position, random_orientation
        class(Abstract_Widom_Method), allocatable :: widom_method
    end type Markov_Chain_Explorer_Wrapper

end module types_markov_chain_explorer_wrapper
