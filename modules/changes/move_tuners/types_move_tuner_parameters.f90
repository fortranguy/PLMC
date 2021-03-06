module types_move_tuner_parameters

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, public :: Concrete_Move_Tuner_Parameters
        integer :: accumulation_period = 0
        real(DP) :: wanted_success_ratio = 0._DP
        real(DP) :: tolerance = 0._DP
    end type Concrete_Move_Tuner_Parameters

end module types_move_tuner_parameters
