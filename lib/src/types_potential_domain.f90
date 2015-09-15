module types_potential_domain

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, public :: Concrete_Potential_Domain
        real(DP) :: min
        real(DP) :: max
        real(DP) :: delta
    end type Concrete_Potential_Domain

end module types_potential_domain
