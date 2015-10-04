module types_ewald_wrapper

use data_constants, only: num_components
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair
use class_ewald_real_particles, only: Abstract_Ewald_Real_Particles

implicit none

private

    type, public :: Ewald_Wrapper
        class(Abstract_Ewald_Real_Pair), allocatable :: real_pair
        class(Abstract_Ewald_Real_Particles), allocatable :: real_particles
    end type Ewald_Wrapper

    type, public :: Inter_Ewald_Wrapper
        class(Abstract_Ewald_Real_Pair), allocatable :: real_pair
    end type Inter_Ewald_Wrapper

    type, public :: Mixture_Ewald_Wrapper
        type(Ewald_Wrapper) :: intras(num_components)
        type(Inter_Ewald_Wrapper) :: inter
    end type Mixture_Ewald_Wrapper

end module types_ewald_wrapper
