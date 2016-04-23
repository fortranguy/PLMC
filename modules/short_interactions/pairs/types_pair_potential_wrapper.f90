module types_pair_potential_wrapper

use classes_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, public :: Pair_Potential_Wrapper
        class(Abstract_Pair_Potential), allocatable :: potential
    end type Pair_Potential_Wrapper

    type, public :: Pair_Potentials_Line
        type(Pair_Potential_Wrapper), allocatable :: line(:)
    end type Pair_Potentials_Line

end module types_pair_potential_wrapper
