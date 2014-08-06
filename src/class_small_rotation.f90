module class_small_rotation

use, intrinsic :: iso_fortran_env, only: DP => REAL64, error_unit
use data_precisions, only: real_zero

implicit none

private

    type, public :: Small_Rotation
    
        real(DP) :: delta
        real(DP) :: delta_save
        real(DP) :: delta_max
        real(DP) :: reject_fix
        
    contains
    
        procedure :: init => Small_Rotation_init
        procedure :: adapt_delta => Small_Rotation_adapt_delta
        procedure :: set_delta => Small_Rotation_set_delta
        procedure :: get_delta => Small_Rotation_get_delta
    
    end type Small_Rotation
    
contains

    pure subroutine Small_Rotation_init(this, delta, delta_max, reject_fix)
    
        class(Small_Rotation), intent(out) :: this
        real(DP), intent(in) :: delta, delta_max
        real(DP), intent(in) :: reject_fix
        
        this%delta = delta
        this%delta_save = this%delta
        this%delta_max = delta_max
        this%reject_fix = reject_fix
        
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
        
        if (reject < this%reject_fix - reject_eps) then
            this%delta = this%delta * more
            if (this%delta > this%delta_max) then
                this%delta = this%delta_max
            end if
        else if (reject > this%reject_fix + reject_eps) then
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
            this%delta = this%delta_save
            write(error_unit, *) "default delta: ", this%delta
        end if
        
        if (this%delta > this%delta_max) then
            write(error_unit, *) type_name, ":   Warning: delta too big."
            this%delta = this%delta_max
            write(error_unit, *) "big delta: ", this%delta
        end if
        
        write(report_unit, *) "Rotation: "
        write(report_unit, *) "    delta = ", this%delta
        write(report_unit, *) "    rejection relative difference = ", &
                                    abs(reject-this%reject_fix)/this%reject_fix
    
    end subroutine Small_Rotation_set_delta

end module class_small_rotation
