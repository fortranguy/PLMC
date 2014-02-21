!> \brief Histogram program of observables

program histogram

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use data_precisions, only : DP
use data_monteCarlo, only : Nstep

implicit none

    ! Physics
    
    integer :: Nobs
    integer :: iStep, iStepIn

    ! Numerical
    
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
    
    allocate(observables(Nobs, Nstep))
    
    write(output_unit, *) "Nstep = ", Nstep
    
    do iStep = 1, Nstep
        read(obs_unit, *) iStepIn, observables(:, iStep)
    end do
    
    open(newunit=histo_unit, recl=4096, file=file(1:length-4)//"_histo.out", action='write')
    
    close(histo_unit)
    
    close(obs_unit)

    deallocate(observables)

end program histogram
