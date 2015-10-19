module types_ewalds_wrapper

use data_constants, only: num_components
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair
use class_ewald_real_component, only: Abstract_Ewald_Real_Component
use class_weighted_structure, only: Abstract_Weighted_Structure

implicit none

private

    type, public :: Ewald_Wrapper
        class(Abstract_Ewald_Real_Pair), allocatable :: real_pair
        class(Abstract_Ewald_Real_Component), allocatable :: real_component
        class(Abstract_Weighted_Structure), allocatable :: weighted_structure
    end type Ewald_Wrapper

    type, public :: Ewald_Wrapper_Macro
        class(Abstract_Ewald_Real_Component), allocatable :: real_component
        class(Abstract_Weighted_Structure), allocatable :: weighted_structure
    end type Ewald_Wrapper_Macro

    type, public :: Ewald_Wrapper_Micro
        class(Abstract_Ewald_Real_Pair), allocatable :: real_pair
    end type Ewald_Wrapper_Micro

    type, public :: Mixture_Ewald_Wrapper
        type(Ewald_Wrapper) :: intras(num_components)
        type(Ewald_Wrapper_Macro) :: inters(num_components)
        type(Ewald_Wrapper_Micro) :: inter_micro
    end type Mixture_Ewald_Wrapper

    type, public :: Ewald_Real_Pair_Wrapper
        class(Abstract_Ewald_Real_Pair), allocatable :: real_pair
    end type Ewald_Real_Pair_Wrapper

    type, public :: Ewald_Real_Pairs_Wrapper
        type(Ewald_Real_Pair_Wrapper), allocatable :: with_components(:)
    end type Ewald_Real_Pairs_Wrapper

    type, public :: Ewalds_Wrapper
        type(Ewald_Real_Pairs_Wrapper), allocatable :: real_pairs(:)
    end type Ewalds_Wrapper

end module types_ewalds_wrapper
