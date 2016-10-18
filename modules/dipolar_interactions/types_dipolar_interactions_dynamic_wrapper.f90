module types_dipolar_interactions_dynamic_wrapper

use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use classes_des_real_component, only: DES_Real_Component_Wrapper
use classes_des_reci_visitor, only: Abstract_DES_Reci_Visitor
use classes_des_self_component, only: DES_Self_Component_Wrapper
use classes_des_surf_mixture, only: Abstract_DES_Surf_Mixture
use classes_dlc_visitor, only: Abstract_DLC_Visitor

implicit none

private

    type, public :: Dipolar_Interactions_Dynamic_Wrapper
        class(Abstract_DES_Convergence_Parameter), allocatable :: alpha
        type(DES_Real_Component_Wrapper), allocatable :: real_components(:, :)
        class(Abstract_DES_Reci_Visitor), allocatable :: reci_visitor
        type(DES_Self_Component_Wrapper), allocatable :: self_components(:)
        class(Abstract_DES_Surf_Mixture), allocatable :: surf_mixture
        class(Abstract_DLC_Visitor), allocatable :: dlc_visitor
    end type Dipolar_Interactions_Dynamic_Wrapper

end module types_dipolar_interactions_dynamic_wrapper
