module class_small_move

use, intrinsic :: iso_fortran_env, only: DP => REAL64, error_unit
use data_precisions, only: real_zero
use data_box, only: num_dimensions
use json_module, only: json_value, json_value_create, to_object, json_value_add

implicit none

private

    type, public :: Small_Move
    
        private
        real(DP), dimension(num_dimensions) :: delta
        real(DP), dimension(num_dimensions) :: delta_save
        real(DP) :: reject_fix
    
    contains
    
        procedure :: init => Small_Move_init
        procedure :: get_delta => Small_Move_get_delta
        procedure :: get_delta_scalar => Small_Move_get_delta_scalar
        procedure :: adapt_delta => Small_Move_adapt_delta
        procedure :: set_delta => Small_Move_set_delta
    
    end type Small_Move
    
contains

    pure subroutine Small_Move_init(this, delta, reject_fix)
    
        class(Small_Move), intent(out) :: this
        real(DP), intent(in) :: delta
        real(DP), intent(in) :: reject_fix
        
        this%delta(:) = delta
        this%delta_save(:) = this%delta(:)
        this%reject_fix = reject_fix
        
    end subroutine Small_Move_init
    
    pure function Small_Move_get_delta(this) result(get_delta)
        class(Small_Move), intent(in) :: this
        real(DP), dimension(num_dimensions) :: get_delta
        
        get_delta(:) = this%delta(:)
    end function Small_Move_get_delta
    
    pure function Small_Move_get_delta_scalar(this) result(get_delta_scalar)
        class(Small_Move), intent(in) :: this
        real(DP) :: get_delta_scalar
        
        get_delta_scalar = sum(this%delta)/size(this%delta)
    end function Small_Move_get_delta_scalar
        
    !> Adapt the displacement delta during thermalisation
    
    pure subroutine Small_Move_adapt_delta(this, Box_size, reject)
    
        class(Small_Move), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size ! warning: average  ?
        real(DP), intent(in) :: reject
        
        real(DP), parameter :: delta_eps = 0.05_DP
        real(DP), parameter :: reject_eps = 0.1_DP * delta_eps
        real(DP), parameter :: more = 1._DP+delta_eps
        real(DP), parameter :: less = 1._DP-delta_eps
        
        if (reject < this%reject_fix - reject_eps) then
            this%delta(:) = this%delta(:) * more
            if (norm2(this%delta) > norm2(Box_size)) then
                this%delta(:) = Box_size(:)
            end if
        else if (reject > this%reject_fix + reject_eps) then
            this%delta(:) = this%delta(:) * less
        end if
    
    end subroutine Small_Move_adapt_delta
    
    subroutine Small_Move_set_delta(this, Box_size, type_name, reject, report_json)
    
        class(Small_Move), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size ! warning: average ?
        character(len=*), intent(in) :: type_name
        real(DP), intent(in) :: reject
        type(json_value), pointer, intent(in) :: report_json

        type(json_value), pointer :: displacement_json

        if (reject < real_zero) then
            write(error_unit, *) type_name, ": Warning: delta adaptation problem."
            this%delta(:) = this%delta_save(:)
            write(error_unit, *) "default delta: ", this%delta(:)
        end if

        if (norm2(this%delta) > norm2(Box_size)) then
            write(error_unit, *) type_name, ": Warning: delta too big."
            this%delta(:) = Box_size(:)
            write(error_unit, *) "big delta: ", this%delta(:)
        end if

        call json_value_create(displacement_json)
        call to_object(displacement_json, "Displacement")
        call json_value_add(report_json, displacement_json)

        call json_value_add(displacement_json, "delta", this%delta)
        call json_value_add(displacement_json, "rejection relative difference", &
                                               abs(reject-this%reject_fix)/this%reject_fix)

        nullify(displacement_json)
    
    end subroutine Small_Move_set_delta

end module class_small_move
