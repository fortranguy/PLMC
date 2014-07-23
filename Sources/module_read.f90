module module_read

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_box, only: num_dimensions

implicit none

private
public jump_snap

contains

    subroutine jump_snap(snap_factors_lcm, snap_factor, num_particles, positions_unit)
    
        integer, intent(in) :: snap_factors_lcm, snap_factor
        integer, intent(in) :: num_particles
        integer, intent(in) :: positions_unit
        
        integer :: snap_factor_jump
        integer :: i_jump
        integer :: i_particle
        real(DP), dimension(num_dimensions) :: dummy_position
        
        snap_factor_jump = snap_factors_lcm / snap_factor

        if (snap_factor_jump > 1) then
            do i_jump = 1, snap_factor_jump - 1
                do i_particle = 1, num_particles
                    read(positions_unit, *) dummy_position(:)
                end do
            end do
        end if
        
    end subroutine jump_snap

end module module_read
