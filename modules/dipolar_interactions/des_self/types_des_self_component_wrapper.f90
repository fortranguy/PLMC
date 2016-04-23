module types_des_self_component_wrapper

use classes_des_self_component, only: Abstract_DES_Self_Component

implicit none

private

    type, public :: DES_Self_Component_Wrapper
        class(Abstract_DES_Self_Component), allocatable :: component
    end type DES_Self_Component_Wrapper

end module types_des_self_component_wrapper
