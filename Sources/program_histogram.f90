!> \brief Histogram program of observables: energy

program histogram

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP
use data_monte_carlo, only: num_equilibrium_steps
use data_histogram, only: hist_dE

implicit none

    integer :: Nobs
    integer :: iStep, iStepIn
    integer :: iObs, iDist
    real(DP), dimension(:, :), allocatable :: limit_values
    integer, dimension(:, :), allocatable :: Ndists
    real(DP), dimension(:), allocatable :: dist_funct_energy
    real(DP) :: x_i
    
    integer :: obs_unit, histo_unit
    character(len=1) :: comment_symbol
    real(DP), dimension(:, :), allocatable :: observables
    
    character(len=4096) :: file
    integer :: length, stat
    logical :: exist
    
    call get_command_argument(1, file, length, stat)
    if (stat /= 0) error stop "error get_command_argument"
    inquire(file=file(1:length), exist=exist)
    if (.not.exist) then
        write(error_unit, *) "missing file: ", file(1:length)
        error stop
    end if
    
    open(newunit=obs_unit, recl=4096, file=file(1:length), status='old', action='read')
    read(obs_unit, *) comment_symbol, Nobs
    write(output_unit, *) "Nobs = ", Nobs
    
    allocate(observables(Nobs, num_equilibrium_steps))
    allocate(limit_values(2, Nobs))
    allocate(Ndists(2, Nobs))
    
    write(output_unit, *) "num_equilibrium_steps = ", num_equilibrium_steps
    
    do iStep = 1, num_equilibrium_steps
        read(obs_unit, *) iStepIn, observables(:, iStep)
    end do
    
    do iObs = 1, Nobs
        limit_values(1, iObs) = minval(observables(iObs, :))
        limit_values(2, iObs) = maxval(observables(iObs, :))
    end do
    Ndists(:, 1) = int(limit_values(:, 1)/hist_dE)
    allocate(dist_funct_energy(Ndists(1, 1): Ndists(2, 1)))
    
    dist_funct_energy(:) = 0._DP
    do iStep = 1, num_equilibrium_steps
        iDist = int(observables(1, iStep)/hist_dE)
        dist_funct_energy(iDist) = dist_funct_energy(iDist) + 1._DP
    end do
    dist_funct_energy(:) = dist_funct_energy(:) / real(num_equilibrium_steps, DP) / hist_dE
    
    open(newunit=histo_unit, recl=4096, file=file(1:length-4)//"_histo.out", action='write')
    do iDist = Ndists(1, 1), Ndists(2, 1)
        x_i = (real(iDist, DP)+0.5_DP) * hist_dE
        write(histo_unit, *) x_i, dist_funct_energy(iDist)
    end do
    
    close(histo_unit)
    close(obs_unit)
    
    deallocate(dist_funct_energy)
    deallocate(Ndists)
    deallocate(limit_values)
    deallocate(observables)

end program histogram
