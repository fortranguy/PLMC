module class_random_orientation

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_orientation, only: gauss

implicit none

private

    type, abstract, public :: Abstract_Random_Orientation
    private
    contains
        procedure(Abstract_Random_Orientation_orientation), deferred, nopass :: orientation
    end type Abstract_Random_Orientation

    abstract interface
        function Abstract_Random_Orientation_orientation() result(orientation)
        import :: DP, num_dimensions
            real(DP) :: orientation(num_dimensions)
        end function Abstract_Random_Orientation_orientation    
    end interface

    type, extends(Abstract_Random_Orientation) :: Null_Random_Orientation
    contains
        procedure, nopass :: orientation => Null_Random_Orientation_orientation
    end type Null_Random_Orientation

    type, extends(Abstract_Random_Orientation) :: Random_Orientation
    contains
        procedure, nopass :: orientation => Random_Orientation_orientation
    end type Random_Orientation

contains

!implementation Null_Random_Orientation
    function Null_Random_Orientation_orientation() result(orientation)
        real(DP) :: orientation(num_dimensions)
        
        orientation = 0._DP
    end function Null_Random_Orientation_orientation
!end implementation Null_Random_Orientation

!implementation Random_Orientation
    !> From SMAC, Algorithm 1.23, p. 43
    function Random_Orientation_orientation() result(orientation)
        real(DP) :: orientation(num_dimensions)

        integer :: i_dimension

        do i_dimension = 1, num_dimensions
            orientation(i_dimension) = gauss()
        end do
        orientation(:) = orientation(:) / norm2(orientation)
    end function Random_Orientation_orientation
!end implementation Random_Orientation

end module class_random_orientation