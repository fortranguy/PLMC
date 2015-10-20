module types_mixture_wrapper

use class_minimum_distance, only: Abstract_Minimum_Distance
use types_component_wrapper, only: Component_Wrapper

implicit none

private

    type, public :: Minimum_Distance_Wrapper
        class(Abstract_Minimum_Distance), allocatable :: min_distance
    end type Minimum_Distance_Wrapper

    type, public :: Minimum_Distances_Wrapper
        type(Minimum_Distance_Wrapper), allocatable :: with_components(:)
    end type Minimum_Distances_Wrapper

    type, public :: Mixture_Wrapper
        type(Component_Wrapper), allocatable :: components(:)
        type(Minimum_Distances_Wrapper), allocatable :: inter_min_distances(:)
        type(Minimum_Distance_Wrapper), allocatable :: wall_min_distances(:)
        logical, allocatable :: components_are_dipolar(:)
    end type Mixture_Wrapper

end module types_mixture_wrapper
