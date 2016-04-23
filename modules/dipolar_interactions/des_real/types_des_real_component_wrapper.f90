module types_des_real_component_wrapper

use classes_des_real_component, only: Abstract_DES_Real_Component

implicit none

private

    type, public :: DES_Real_Component_Wrapper
        class(Abstract_DES_Real_Component), allocatable :: component
    end type DES_Real_Component_Wrapper

end module types_des_real_component_wrapper
