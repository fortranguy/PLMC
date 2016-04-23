module types_des_real_pair_wrapper

use classes_des_real_pair, only: Abstract_DES_Real_Pair

implicit none

private

    type, public :: DES_Real_Pair_Wrapper
        class(Abstract_DES_Real_Pair), allocatable :: potential
    end type DES_Real_Pair_Wrapper

    type, public :: DES_Real_Pairs_Line
        type(DES_Real_Pair_Wrapper), allocatable :: line(:)
    end type DES_Real_Pairs_Line

end module types_des_real_pair_wrapper
