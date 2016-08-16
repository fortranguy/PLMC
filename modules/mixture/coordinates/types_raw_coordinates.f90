module types_raw_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, public :: Concrete_Raw_Coordinates
        real(DP), dimension(:, :), allocatable :: positions, orientations
    end type Concrete_Raw_Coordinates

end module types_raw_coordinates
