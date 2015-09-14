module types_potential_parameters

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, abstract, public :: Abstract_Potential_Parameters
    end type Abstract_Potential_Parameters

    type, extends(Abstract_Potential_Parameters), public :: Null_Potential_Paramters
    end type Null_Potential_Paramters

    type, extends(Abstract_Potential_Parameters), public :: Lennard_Jones_Parameters
        real(DP) :: epsilon
        real(DP) :: sigma
    end type Lennard_Jones_Parameters

end module types_potential_parameters

