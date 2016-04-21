module types_mixture_wrapper

use types_component_wrapper, only: Component_Wrapper
use types_min_distances_wrapper, only: Min_Distance_Wrapper, Min_Distances_Wrapper
use classes_mixture_total_moment, only: Abstract_Mixture_Total_Moment

implicit none

private

    type, public :: Mixture_Wrapper
        type(Component_Wrapper), allocatable :: components(:)
        type(Min_Distances_Wrapper), allocatable :: components_min_distances(:)
        type(Min_Distance_Wrapper), allocatable :: wall_min_distances(:)
        class(Abstract_Mixture_Total_Moment), allocatable :: total_moment
    end type Mixture_Wrapper

end module types_mixture_wrapper
