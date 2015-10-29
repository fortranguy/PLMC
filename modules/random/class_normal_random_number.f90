module class_normal_random_number

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    real(DP), parameter :: sigma3d = 1._DP / sqrt(3._DP)

    type, public :: Normal_Random_Number
    private
        logical :: switch = .false.
        real(DP) :: saved_result
    contains
        procedure :: get => Normal_Random_Number_get
    end type Normal_Random_Number

contains

    !> From SMAC, Algorithm 1.19, p.39
    function Normal_Random_Number_get(this) result(gauss)
        class(Normal_Random_Number), intent(inout) :: this
        real(DP) :: gauss

        real(DP), dimension(2) :: xy, rand2d
        real(DP) :: norm_sqr, factor

        if (.not.this%switch) then
            do
                call random_number(rand2d)
                xy = 2._DP*rand2d - 1._DP
                norm_sqr = sum(xy**2)
                if (0._DP < norm_sqr .and. norm_sqr <= 1._DP) exit
            end do
            factor = sqrt(-2._DP * log(norm_sqr) / norm_sqr)
            gauss = sigma3d * factor * xy(1)
            this%saved_result = sigma3d * factor * xy(2)
            this%switch = .true.
        else
            gauss = this%saved_result
            this%switch = .false.
        end if
    end function Normal_Random_Number_get

end module class_normal_random_number
