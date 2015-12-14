module types_long_interactions_wrapper

use class_ewald_convergence_parameter, only: Abstract_Ewald_Convergence_Parameter
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair
use class_ewald_real_component, only: Abstract_Ewald_Real_Component
use class_ewald_real_visitor, only: Abstract_Ewald_Real_Visitor
use class_ewald_reci_weight, only: Abstract_Ewald_Reci_Weight
use class_ewald_reci_structure_, only: Abstract_Ewald_Reci_Structure_
use class_ewald_reci_delta_visitor, only: Abstract_Ewald_Reci_Delta_Visitor
use class_ewald_self_component, only: Abstract_Ewald_Self_Component

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

    type, public :: Ewald_Reci_Component_Wrapper
        class(Abstract_Ewald_Reci_Structure_), allocatable :: reci_structure
        class(Abstract_Ewald_Reci_Delta_Visitor), allocatable :: reci_delta_visitor
    end type Ewald_Reci_Component_Wrapper

    type, public :: Ewald_Self_Component_Wrapper
        class(Abstract_Ewald_Self_Component), allocatable :: self
    end type Ewald_Self_Component_Wrapper

    !Dipolar interactions only
    type, public :: Long_Interactions_Wrapper
        class(Abstract_Ewald_Convergence_Parameter), allocatable :: alpha
        class(Abstract_Ewald_Real_Visitor), allocatable :: real_visitor
        type(Ewald_Real_Component_Wrapper), allocatable :: real_components(:, :)
        type(Ewald_Real_Pairs_Wrapper), allocatable :: real_pairs(:)
        class(Abstract_Ewald_Reci_Weight), allocatable :: reci_weight
        type(Ewald_Reci_Component_Wrapper), allocatable :: reci_components(:)
        type(Ewald_Self_Component_Wrapper), allocatable :: self_components(:)
    end type Long_Interactions_Wrapper

end module types_long_interactions_wrapper
