module module_read

use data_precisions, only: DP
use data_box, only: Ndim, Height
use data_distribution, only: z_min_ratio, z_max_ratio

implicit none

private
public jump_snap, test_particle_inside

contains

    subroutine jump_snap(snap_factors_lcm, snap_factor, Ncol, positions_unit)
    
        integer, intent(in) :: snap_factors_lcm, snap_factor
        integer, intent(in) :: Ncol
        integer, intent(in) :: positions_unit
        
        integer :: snap_factor_jump
        integer :: iJump
        integer :: iCol
        real(DP), dimension(Ndim) :: dummy_position
        
        snap_factor_jump = snap_factors_lcm / snap_factor

        if (snap_factor_jump > 1) then
            do iJump = 1, snap_factor_jump - 1
                do iCol = 1, Ncol
                    read(positions_unit, *) dummy_position(:)
                end do
            end do
        end if
        
    end subroutine jump_snap
    
    subroutine test_particle_inside(Ncol, positions, paricles_inside, Ncol_step)
    
        integer, intent(in) :: Ncol
        real(DP), dimension(:, :), intent(in) :: positions
        logical, dimension(:), intent(out) :: paricles_inside
        integer, intent(out) :: Ncol_step
        
        integer :: iCol
        
        Ncol_step = 0
        do iCol = 1, Ncol        
            if (z_min_ratio * Height < positions(3, iCol) .and. &
                positions(3, iCol) < z_max_ratio * Height) then
                paricles_inside(iCol) = .true.
                Ncol_step = Ncol_step + 1
            else
                paricles_inside(iCol) = .false.
            end if
        end do
        
    end subroutine test_particle_inside

end module module_read

!> \brief Calculate and print the ``mixing distribution'' function

program mix_distribution

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP, real_zero
use data_box, only: Ndim, Lsize, Height
use data_monteCarlo, only: Nstep
use data_distribution, only: snap, dist_dr, z_min_ratio, z_max_ratio
use module_maths, only: lcm
use module_physics_micro, only: sphere_volume, dist_PBC
use module_arguments, only: argument_to_file
use module_read, only: jump_snap, test_particle_inside
!$ use omp_lib

implicit none
   
    real(DP) :: Volume_inside

    character(len=5) :: type1_name, type2_name
    integer :: type1_Ncol, type2_Ncol
    integer :: type1_snap_factor, type2_snap_factor
    integer :: type1_positions_unit, type2_positions_unit
    real(DP) :: type1_density_inside, type2_density_inside
    logical, dimension(:), allocatable :: type1_paricles_inside, type2_paricles_inside
    real(DP) :: type1_Ncol_inside, type2_Ncol_inside
    integer :: type1_Ncol_sum, type2_Ncol_sum
    integer :: type1_Ncol_step, type2_Ncol_step
    integer :: snap_factors_lcm
    
    integer :: report_unit, distrib_unit
    integer :: Ndist
    integer :: iStep
    integer :: type1_iCol, type2_iCol
    integer :: iDist
    real(DP) :: distance_12
    real(DP) :: distance_i, distance_minus, distance_plus
    real(DP), parameter :: distance_max = 3._DP
    real(DP), dimension(:), allocatable :: dist_step
    real(DP), dimension(:), allocatable :: dist_function
    real(DP), dimension(:, :), allocatable :: type1_positions, type2_positions
    
    character(len=4096) :: file
    integer :: length
    
    real(DP) :: initial_time, final_time
    
    if (.not.snap) stop "Snap désactivé."
    
    Volume_inside = product(Lsize(1:2)) * (z_max_ratio - z_min_ratio) * (Height-1._DP)
    
    Ndist = int(distance_max / dist_dr)
    allocate(dist_step(Ndist))
    allocate(dist_function(Ndist))
    
    call argument_to_file(1, file, length)
    open(newunit=type1_positions_unit, recl=4096, file=file(1:length), status='old', action='read')    
    read(type1_positions_unit, *) type1_name, type1_Ncol, type1_snap_factor
    write(output_unit, *) "type 1: ", type1_name, type1_Ncol, type1_snap_factor
    allocate(type1_positions(Ndim, type1_Ncol))
    allocate(type1_paricles_inside(type1_Ncol))
    
    call argument_to_file(2, file, length)
    open(newunit=type2_positions_unit, recl=4096, file=file(1:length), status='old', action='read')    
    read(type2_positions_unit, *) type2_name, type2_Ncol, type2_snap_factor
    write(output_unit, *) "type 1: ", type2_name, type2_Ncol, type2_snap_factor
    allocate(type2_positions(Ndim, type2_Ncol))
    allocate(type2_paricles_inside(type2_Ncol))
    
    snap_factors_lcm = lcm(type1_snap_factor, type2_snap_factor)
    
    type1_Ncol_sum = 0
    type2_Ncol_sum = 0
    
    write(output_unit, *) "Start !"
    call cpu_time(initial_time)
    do iStep = 1, Nstep / snap_factors_lcm
        
        call jump_snap(snap_factors_lcm, type1_snap_factor, type1_Ncol, type1_positions_unit)
        do type1_iCol = 1, type1_Ncol
            read(type1_positions_unit, *) type1_positions(:, type1_iCol)
        end do

        call jump_snap(snap_factors_lcm, type2_snap_factor, type2_Ncol, type2_positions_unit)
        do type2_iCol = 1, type2_Ncol
            read(type2_positions_unit, *) type2_positions(:, type2_iCol)
        end do
        
        call test_particle_inside(type1_Ncol, type1_positions, type1_paricles_inside, type1_Ncol_step)
        type1_Ncol_sum = type1_Ncol_sum + type1_Ncol_step
        call test_particle_inside(type2_Ncol, type2_positions, type2_paricles_inside, type2_Ncol_step)
        type2_Ncol_sum = type2_Ncol_sum + type2_Ncol_step
        
        dist_step(:) = 0._DP
        
        do type1_iCol = 1, type1_Ncol
            if (type1_paricles_inside(type1_iCol))  then
            
                do type2_iCol = 1, type2_Ncol                
                    distance_12 = dist_PBC(type1_positions(:, type1_iCol), type2_positions(:, type2_iCol))
                    if (distance_12 <= distance_max) then
                        iDist = int(distance_12 / dist_dr)
                        dist_step(iDist) = dist_step(iDist) + 1._DP
                    end if
                end do
                
            end if
        end do
        
        dist_step(:) = dist_step(:) / real(type1_Ncol_step, DP)
        dist_function(:) = dist_function(:) + dist_step(:)
    
    end do
    call cpu_time(final_time)
    write(output_unit, *) "Finish !"
    
    deallocate(type2_paricles_inside)
    deallocate(type2_positions) 
    close(type2_positions_unit)   
    
    deallocate(type1_paricles_inside)
    deallocate(type1_positions) 
    close(type1_positions_unit)
    
    open(newunit=report_unit, file=type1_name//"-"//type2_name//"_mix_distribution_report.txt", &
         action="write")
    
    write(report_unit, *) "Volume_inside =", Volume_inside
    write(report_unit, *) "z-domain = ", z_min_ratio * Height, z_max_ratio * Height
    
    type1_Ncol_inside = real(type1_Ncol_sum, DP) / real(Nstep/snap_factors_lcm, DP)
    type1_density_inside = type1_Ncol_inside / Volume_inside
    write(report_unit, *) type1_name, " inside density: ", type1_density_inside
    type2_Ncol_inside = real(type2_Ncol_sum, DP) / real(Nstep/snap_factors_lcm, DP)
    type2_density_inside = type2_Ncol_inside / Volume_inside
    write(report_unit, *) type2_name, " inside density: ", type2_density_inside
    write(report_unit, *) "Duration =", (final_time - initial_time) / 60._DP, "min"
    
    close(report_unit)
    
    open(newunit=distrib_unit, file=type1_name//"-"//type2_name//"_mix_distribution.out", &
         action="write")
    
        dist_function(:) = dist_function(:) / real(Nstep/snap_factors_lcm, DP)
    
        do iDist = 1, Ndist
        
            distance_i = (real(iDist, DP) + 0.5_DP) * dist_dr
            distance_minus = real(iDist, DP) * dist_dr
            distance_plus = real(iDist + 1, DP) * dist_dr
            
            dist_function(iDist) = dist_function(iDist) / type2_density_inside / &
                                   (sphere_volume(distance_plus) - sphere_volume(distance_minus))
            write(distrib_unit, *) distance_i, dist_function(iDist)
            
        end do
        
    close(distrib_unit)
    
    deallocate(dist_function)
    deallocate(dist_step)

end program mix_distribution
