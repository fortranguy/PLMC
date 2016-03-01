module types_line_observables

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, public :: Concrete_Line_Observables
        real(DP), allocatable :: with_components(:)
    end type Concrete_Line_Observables

end module types_line_observables
