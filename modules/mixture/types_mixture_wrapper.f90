module types_mixture_wrapper

use classes_minimum_distance, only: Abstract_Minimum_Distance
use types_component_wrapper, only: Component_Wrapper
use classes_mixture_total_moment, only: Abstract_Mixture_Total_Moment

implicit none

private

    type, public :: Minimum_Distance_Wrapper
        class(Abstract_Minimum_Distance), allocatable :: distance
    end type Minimum_Distance_Wrapper

    type, public :: Minimum_Distances_Wrapper
        type(Minimum_Distance_Wrapper), allocatable :: line(:)
    end type Minimum_Distances_Wrapper

    type, public :: Mixture_Wrapper
        type(Component_Wrapper), allocatable :: components(:)
        type(Minimum_Distances_Wrapper), allocatable :: components_min_distances(:)
        type(Minimum_Distance_Wrapper), allocatable :: wall_min_distances(:)
        class(Abstract_Mixture_Total_Moment), allocatable :: total_moment
    end type Mixture_Wrapper

end module types_mixture_wrapper
