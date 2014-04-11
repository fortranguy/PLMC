module class_small_move

use, intrinsic :: iso_fortran_env, only: error_unit
use data_precisions, only: DP, real_zero
use data_box, only: Ndim

implicit none

private

    type, public :: Small_move
    
        real(DP), dimension(Ndim) :: delta
        real(DP), dimension(Ndim) :: deltaSave
        real(DP) :: rejectFix
    
    contains
    
        procedure :: init => Small_move_init
        procedure :: adapt_delta => Small_move_adapt_delta
        procedure :: set_delta => Small_move_set_delta
    
    end type Small_move
    
contains

    pure subroutine Small_move_init(this, delta, rejectFix)
        class(Small_move), intent(out) :: this
        real(DP), dimension(Ndim), intent(in) :: delta
        real(DP), intent(in) :: rejectFix
        
        this%delta = delta
        this%deltaSave = this%delta
        this%rejectFix = rejectFix
    end subroutine Small_move_init
        
    !> Adapt the displacement delta during thermalisation
    
    pure subroutine Small_move_adapt_delta(this, Box_size, reject)
    
        class(Small_move), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size ! warning: average  ?
        real(DP), intent(in) :: reject
        
        real(DP), parameter :: delta_eps = 0.05_DP
        real(DP), parameter :: move_reject_eps = 0.1_DP * delta_eps
        real(DP), parameter :: more = 1._DP+delta_eps
        real(DP), parameter :: less = 1._DP-delta_eps
        
        if (reject < this%rejectFix - move_reject_eps) then
            this%delta(:) = this%delta(:) * more
            if (norm2(this%delta) > norm2(Box_size)) then
                this%delta(:) = Box_size(:)
            end if
        else if (reject > this%rejectFix + move_reject_eps) then
            this%delta(:) = this%delta(:) * less
        end if
    
    end subroutine Small_move_adapt_delta
    
    subroutine Small_move_set_delta(this, type_name, Box_size, reject, report_unit)
    
        class(Small_move), intent(inout) :: this
        character(len=*), intent(in) :: type_name
        real(DP), dimension(:), intent(in) :: Box_size ! warning: average ?
        real(DP), intent(in) :: reject
        integer, intent(in) :: report_unit

        if (reject < real_zero) then
            write(error_unit, *) type_name, ":    Warning: delta adaptation problem."
            this%delta(:) = this%deltaSave(:)
            write(error_unit, *) "default delta: ", this%delta(:)
        end if

        if (norm2(this%delta) > norm2(Box_size)) then
            write(error_unit, *) type_name, ":   Warning: delta too big."
            this%delta(:) = Box_size(:)
            write(error_unit, *) "big delta: ", this%delta(:)
        end if

        write(report_unit, *) "Displacement: "
        write(report_unit, *) "    delta(:) = ", this%delta(:)
        write(report_unit, *) "    rejection relative difference = ", &
                                    abs(reject-this%rejectFix)/this%rejectFix
    
    end subroutine Small_move_set_delta

end module class_small_move
