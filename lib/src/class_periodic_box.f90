module class_periodic_box

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use data_precisions, only: real_zero
use procedures_errors, only: warning_continue, error_exit
use procedures_checks, only: check_3d_array, check_positive

implicit none

private

    type, abstract, public :: Abstract_Periodic_Box
    private
        real(DP) :: real_size(num_dimensions)
    contains
        procedure :: set_real_size => Abstract_Periodic_Box_set_real_size
        procedure, nopass, private :: check => Abstract_Periodic_Box_check
        procedure :: get_real_size => Abstract_Periodic_Box_get_real_size
        procedure :: distance => Abstract_Periodic_Box_distance
        procedure(Abstract_Periodic_Box_vector), deferred :: vector
    end type Abstract_Periodic_Box

    abstract interface
    
        pure function Abstract_Periodic_Box_vector(this, position_1, position_2) result(vector)
        import :: DP, num_dimensions, Abstract_Periodic_Box
            class(Abstract_Periodic_Box), intent(in) :: this
            real(DP), intent(in) :: position_1(:), position_2(:)
            real(DP) :: vector(num_dimensions)
        end function Abstract_Periodic_Box_vector
        
    end interface

    type, extends(Abstract_Periodic_Box), public :: XYZ_Periodic_Box
    contains
        procedure :: vector => XYZ_Periodic_Box_vector
    end type XYZ_Periodic_Box

    type, extends(Abstract_Periodic_Box), public :: XY_Periodic_Box
    contains
        procedure :: vector => XY_Periodic_Box_vector
    end type XY_Periodic_Box
    
contains

!implementation Abstract_Periodic_Box

    subroutine Abstract_Periodic_Box_set_real_size(this, real_size)
        class(Abstract_Periodic_Box), intent(out) :: this
        real(DP), intent(in) :: real_size(:)

        call this%check(real_size)
        this%real_size = real_size
    end subroutine Abstract_Periodic_Box_set_real_size

    subroutine Abstract_Periodic_Box_check(real_size)
        real(DP), intent(in) :: real_size(:)
        
        call check_3d_array("Abstract_Periodic_Box", "real_size", real_size)
        call check_positive("Abstract_Periodic_Box", "real_size", real_size)
        if (abs(real_size(1)-real_size(2)) > real_zero) then
            call warning_continue("Abstract_Periodic_Box: "//&
                "real_size(1) and real_size(2) are not equal.")
        end if
    end subroutine Abstract_Periodic_Box_check

    pure function Abstract_Periodic_Box_get_real_size(this) result(real_size)
        class(Abstract_Periodic_Box), intent(in) :: this
        real(DP) :: real_size(num_dimensions)

        real_size = this%real_size
    end function Abstract_Periodic_Box_get_real_size

    pure function Abstract_Periodic_Box_distance(this, position_1, position_2) result(distance)
        class(Abstract_Periodic_Box), intent(in) :: this
        real(DP), intent(in) :: position_1(:), position_2(:)
        real(DP) :: distance

        distance = norm2(this%vector(position_1, position_2))
    end function Abstract_Periodic_Box_distance
    
!end implementation Abstract_Periodic_Box

!implementation XYZ_Periodic_Box

    !> from SMAC, algorithm 2.5 & 2.6, p.91
    pure function XYZ_Periodic_Box_vector(this, position_1, position_2) result(vector)
        class(XYZ_Periodic_Box), intent(in) :: this
        real(DP), intent(in) :: position_1(:), position_2(:)
        real(DP) :: vector(num_dimensions)

        vector = modulo(position_2 - position_1, this%real_size)

        where(vector > this%real_size/2._DP)
            vector = vector - this%real_size
        end where
    end function XYZ_Periodic_Box_vector

!end implementation XYZ_Periodic_Box

!implementation XY_Periodic_Box

    pure function XY_Periodic_Box_vector(this, position_1, position_2) result(vector)
        class(XY_Periodic_Box), intent(in) :: this
        real(DP), intent(in) :: position_1(:), position_2(:)
        real(DP) :: vector(num_dimensions)

        vector(1:2) = modulo(position_2(1:2) - position_1(1:2), this%real_size(1:2))

        where(vector(1:2) > this%real_size(1:2)/2._DP)
            vector(1:2) = vector(1:2) - this%real_size(1:2)
        end where

        vector(3) = position_2(3) - position_1(3)
    end function XY_Periodic_Box_vector

!end implementation XY_Periodic_Box

end module class_periodic_box
