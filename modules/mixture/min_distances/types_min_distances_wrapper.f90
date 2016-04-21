module types_min_distances_wrapper

use classes_min_distance, only: Abstract_Min_Distance

implicit none

private

    type, public :: Min_Distance_Wrapper
        class(Abstract_Min_Distance), allocatable :: distance
    end type Min_Distance_Wrapper

    type, public :: Min_Distances_Wrapper
        type(Min_Distance_Wrapper), allocatable :: line(:)
    end type Min_Distances_Wrapper

end module types_min_distances_wrapper
