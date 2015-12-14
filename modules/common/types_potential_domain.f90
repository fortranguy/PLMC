module types_potential_domain

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, abstract, public :: Abstract_Potential_Domain
        real(DP) :: min = 0._DP
        real(DP) :: delta = 0._DP
    end type Abstract_Potential_Domain

    type, extends(Abstract_Potential_Domain), public :: Short_Potential_Domain
        real(DP) :: max = 0._DP
    end type Short_Potential_Domain

    type, extends(Abstract_Potential_Domain), public :: Long_Potential_Domain
        real(DP) :: max_over_box = 0._DP ! volume dependency
    end type Long_Potential_Domain

end module types_potential_domain
