module types_dipolar_interactions_static_wrapper

use classes_box_volume_memento, only: Abstract_Box_Volume_Memento
use classes_des_real_pair, only: Abstract_DES_Real_Pair
use classes_des_reci_weight, only: Abstract_DES_Reci_Weight
use classes_des_reci_structure, only: Abstract_DES_Reci_Structure
use classes_dlc_weight, only: Abstract_DLC_Weight
use classes_dlc_structures, only: Abstract_DLC_Structures

implicit none

private

    type, public :: Dipolar_Interactions_Static_Wrapper
        class(Abstract_Box_Volume_Memento), allocatable :: box_volume_memento_real
        class(Abstract_DES_Real_Pair), allocatable :: real_pair
        class(Abstract_Box_Volume_Memento), allocatable :: box_volume_memento_reci
        class(Abstract_DES_Reci_Weight), allocatable :: reci_weight
        class(Abstract_DES_Reci_Structure), allocatable :: reci_structure
        class(Abstract_DLC_Weight), allocatable :: dlc_weight
        class(Abstract_DLC_Structures), allocatable :: dlc_structures
    end type Dipolar_Interactions_Static_Wrapper

end module types_dipolar_interactions_static_wrapper
