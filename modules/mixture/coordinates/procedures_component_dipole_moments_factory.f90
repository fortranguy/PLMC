module procedures_component_dipole_moments_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments, &
    Concrete_Component_Dipole_Moments, Null_Component_Dipole_Moments

implicit none

private
public :: create, destroy

contains

    subroutine create(dipole_moments, orientations, is_dipolar, generating_data, prefix)
        class(Abstract_Component_Dipole_Moments), allocatable, intent(out) :: dipole_moments
        class(Abstract_Component_Coordinates), intent(in) :: orientations
        logical, intent(in) :: is_dipolar
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: moment_norm

        if (is_dipolar) then
            data_field = prefix//"moment norm"
            call generating_data%get(data_field, moment_norm, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Component_Dipole_Moments :: dipole_moments)
        else
            moment_norm = 0._DP
            allocate(Null_Component_Dipole_Moments :: dipole_moments)
        end if
        call dipole_moments%construct(moment_norm, orientations)
    end subroutine create

    subroutine destroy(dipole_moments)
        class(Abstract_Component_Dipole_Moments), allocatable, intent(inout) :: dipole_moments

        if (allocated(dipole_moments)) then
            call dipole_moments%destroy()
            deallocate(dipole_moments)
        end if
    end subroutine destroy

end module procedures_component_dipole_moments_factory
