module types_dipolar_interactions_wrapper

use class_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use class_des_real_pair, only: Abstract_DES_Real_Pair
use class_des_real_component, only: Abstract_DES_Real_Component
use class_des_real_visitor, only: Abstract_DES_Real_Visitor
use class_des_reci_weight, only: Abstract_DES_Reci_Weight
use class_des_reci_structure, only: Abstract_DES_Reci_Structure
use class_des_reci_visitor, only: Abstract_DES_Reci_Visitor
use class_des_self_component, only: Abstract_DES_Self_Component
use class_des_surf_mixture, only: Abstract_DES_Surf_Mixture
use class_dlc_weight, only: Abstract_DLC_Weight
use class_dlc_structures, only: Abstract_DLC_Structures
use class_dlc_visitor, only: Abstract_DLC_Visitor

implicit none

private

    type, public :: DES_Real_Pair_Wrapper
        class(Abstract_DES_Real_Pair), allocatable :: potential
    end type DES_Real_Pair_Wrapper

    type, public :: DES_Real_Pairs_Wrapper
        type(DES_Real_Pair_Wrapper), allocatable :: line(:)
    end type DES_Real_Pairs_Wrapper

    type, public :: DES_Real_Component_Wrapper
        class(Abstract_DES_Real_Component), allocatable :: component
    end type DES_Real_Component_Wrapper

    type, public :: DES_Self_Component_Wrapper
        class(Abstract_DES_Self_Component), allocatable :: component
    end type DES_Self_Component_Wrapper

    type, public :: Dipolar_Interactions_Wrapper
        class(Abstract_DES_Convergence_Parameter), allocatable :: alpha
        class(Abstract_DES_Real_Visitor), allocatable :: real_visitor
        type(DES_Real_Component_Wrapper), allocatable :: real_components(:, :)
        type(DES_Real_Pairs_Wrapper), allocatable :: real_pairs(:)
        class(Abstract_DES_Reci_Weight), allocatable :: reci_weight
        class(Abstract_DES_Reci_Structure), allocatable :: reci_structure
        class(Abstract_DES_Reci_Visitor), allocatable :: reci_visitor
        type(DES_Self_Component_Wrapper), allocatable :: self_components(:)
        class(Abstract_DES_Surf_Mixture), allocatable :: surf_mixture
        class(Abstract_DLC_Weight), allocatable :: dlc_weight
        class(Abstract_DLC_Structures), allocatable :: dlc_structures
        class(Abstract_DLC_Visitor), allocatable :: dlc_visitor
    end type Dipolar_Interactions_Wrapper

end module types_dipolar_interactions_wrapper
