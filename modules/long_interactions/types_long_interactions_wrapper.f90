module types_long_interactions_wrapper

use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair
use class_ewald_real_component, only: Abstract_Ewald_Real_Component
use class_weighted_structure, only: Abstract_Weighted_Structure
use class_ewald_real_visitor, only: Abstract_Ewald_Real_Visitor

implicit none

private

    type, public :: Ewald_Real_Pair_Wrapper
        class(Abstract_Ewald_Real_Pair), allocatable :: real_pair
    end type Ewald_Real_Pair_Wrapper

    type, public :: Ewald_Real_Pairs_Wrapper
        type(Ewald_Real_Pair_Wrapper), allocatable :: with_components(:)
    end type Ewald_Real_Pairs_Wrapper

    type, public :: Ewald_Real_Component_Wrapper
        class(Abstract_Ewald_Real_Component), allocatable :: real_component
    end type Ewald_Real_Component_Wrapper

    !Dipolar interactions only
    type, public :: Long_Interactions_Wrapper
        class(Abstract_Ewald_Real_Visitor), allocatable :: real_visitor
        type(Ewald_Real_Component_Wrapper), allocatable :: real_components(:, :)
        type(Ewald_Real_Pairs_Wrapper), allocatable :: real_pairs(:)
    end type Long_Interactions_Wrapper

end module types_long_interactions_wrapper