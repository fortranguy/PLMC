module types_dipolar_interactions_wrapper

use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use classes_des_real_pair, only: Abstract_DES_Real_Pair
use types_des_real_component_wrapper, only: DES_Real_Component_Wrapper
use classes_des_reci_weight, only: Abstract_DES_Reci_Weight
use classes_des_reci_structure, only: Abstract_DES_Reci_Structure
use classes_des_reci_visitor, only: Abstract_DES_Reci_Visitor
use types_des_self_component_wrapper, only: DES_Self_Component_Wrapper
use classes_des_surf_mixture, only: Abstract_DES_Surf_Mixture
use classes_dlc_weight, only: Abstract_DLC_Weight
use classes_dlc_structures, only: Abstract_DLC_Structures
use classes_dlc_visitor, only: Abstract_DLC_Visitor

implicit none

private

    type, public :: Dipolar_Interactions_Wrapper
        class(Abstract_DES_Convergence_Parameter), allocatable :: alpha
        class(Abstract_DES_Real_Pair), allocatable :: real_pair
        type(DES_Real_Component_Wrapper), allocatable :: real_components(:, :)
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
