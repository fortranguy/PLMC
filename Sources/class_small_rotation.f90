module class_small_rotation

use, intrinsic :: iso_fortran_env, only: error_unit
use data_precisions, only: DP, real_zero

implicit none

private

    type, public :: Small_Rotation
    
        real(DP) :: delta
        real(DP) :: deltaSave
        real(DP) :: deltaMax
        real(DP) :: rejectFix
        
    contains
    
        procedure :: init => Small_Rotation_init
        procedure :: adapt_delta => Small_Rotation_adapt_delta
        procedure :: set_delta => Small_Rotation_set_delta
        procedure :: get_delta => Small_Rotation_get_delta
    
    end type Small_Rotation
    
contains

    pure subroutine Small_Rotation_init(this, delta, deltaMax, rejectFix)
        class(Small_Rotation), intent(out) :: this
        real(DP), intent(in) :: delta, deltaMax
        real(DP), intent(in) :: rejectFix
        
        this%delta = delta
        this%deltaSave = this%delta
        this%deltaMax = deltaMax
        this%rejectFix = rejectFix    
    end subroutine Small_Rotation_init

    pure function Small_Rotation_get_delta(this) result(get_delta)
        class(Small_Rotation), intent(in) :: this
        real(DP) :: get_delta
        
        get_delta = this%delta
    end function Small_Rotation_get_delta
    
    !> Adaptation of delta during the thermalisation
    
    pure subroutine Small_Rotation_adapt_delta(this, reject)
    
        class(Small_Rotation), intent(inout) :: this
        real(DP), intent(in) :: reject
        
        real(DP), parameter :: delta_eps = 0.05_DP
        real(DP), parameter :: reject_eps = 0.1_DP * delta_eps
        real(DP), parameter :: more = 1._DP+delta_eps
        real(DP), parameter :: less = 1._DP-delta_eps
        
        if (reject < this%rejectFix - reject_eps) then
            this%delta = this%delta * more
            if (this%delta > this%deltaMax) then
                this%delta = this%deltaMax
            end if
        else if (reject > this%rejectFix + reject_eps) then
            this%delta = this%delta * less
        end if
    
    end subroutine Small_Rotation_adapt_delta
    
    subroutine Small_Rotation_set_delta(this, type_name, reject, report_unit)
    
        class(Small_Rotation), intent(inout) :: this
        character(len=*), intent(in) :: type_name
        real(DP), intent(in) :: reject
        integer, intent(in) :: report_unit
        
        if (reject < real_zero) then
            write(error_unit, *) type_name, ":    Warning: delta adaptation problem."
            this%delta = this%deltaSave
            write(error_unit, *) "default delta: ", this%delta
        end if
        
        if (this%delta > this%deltaMax) then
            write(error_unit, *) type_name, ":   Warning: delta too big."
            this%delta = this%deltaMax
            write(error_unit, *) "big delta: ", this%delta
        end if
        
        write(report_unit, *) "Rotation: "
        write(report_unit, *) "    delta = ", this%delta
        write(report_unit, *) "    rejection relative difference = ", &
                                    abs(reject-this%rejectFix)/this%rejectFix
    
    end subroutine Small_Rotation_set_delta

end module class_small_rotation
