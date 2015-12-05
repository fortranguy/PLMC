module types_long_interactions_wrapper

use class_ewald_convergence_parameter, only: Abstract_Ewald_Convergence_Parameter
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair
use class_ewald_real_component, only: Abstract_Ewald_Real_Component
use class_ewald_reci_structure, only: Abstract_Ewald_Reci_Structure
use class_ewald_real_visitor, only: Abstract_Ewald_Real_Visitor
use class_ewald_reci_weight, only: Abstract_Ewald_Reci_Weight
use class_ewald_reci_structure, only: Abstract_Ewald_Reci_Structure
use class_ewald_self, only: Abstract_Ewald_Self

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

    type, public :: Ewald_Reci_Structure_Wrapper
        class(Abstract_Ewald_Reci_Structure), allocatable :: reci_structure
    end type Ewald_Reci_Structure_Wrapper

    type, public :: Ewald_Self_Wrapper
        class(Abstract_Ewald_Self), allocatable :: self
    end type Ewald_Self_Wrapper

    !Dipolar interactions only
    type, public :: Long_Interactions_Wrapper
        class(Abstract_Ewald_Convergence_Parameter), allocatable :: alpha
        class(Abstract_Ewald_Real_Visitor), allocatable :: real_visitor
        type(Ewald_Real_Component_Wrapper), allocatable :: real_components(:, :)
        type(Ewald_Real_Pairs_Wrapper), allocatable :: real_pairs(:)
        class(Abstract_Ewald_Reci_Weight), allocatable :: reci_weight
        type(Ewald_Reci_Structure_Wrapper), allocatable :: reci_structures(:)
        type(Ewald_Self_Wrapper), allocatable :: self(:)
    end type Long_Interactions_Wrapper

end module types_long_interactions_wrapper
