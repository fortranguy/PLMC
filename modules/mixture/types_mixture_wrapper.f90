module types_mixture_wrapper

use types_component_wrapper, only: Component_Wrapper
use classes_average_num_particles, only: Abstract_Average_Num_Particles
use classes_min_distance, only: Min_Distance_Wrapper, Min_Distance_Line
use classes_mixture_total_moment, only: Abstract_Mixture_Total_Moment

implicit none

private

    type, public :: Mixture_Wrapper
        type(Component_Wrapper), allocatable :: components(:, :)
        class(Abstract_Average_Num_Particles), allocatable :: average_nums_particles(:, :)
        type(Min_Distance_Line), allocatable :: components_min_distances(:)
        type(Min_Distance_Wrapper), allocatable :: wall_min_distances(:)
        class(Abstract_Mixture_Total_Moment), allocatable :: total_moments(:)
    end type Mixture_Wrapper

end module types_mixture_wrapper
