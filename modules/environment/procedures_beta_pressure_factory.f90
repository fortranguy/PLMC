module procedures_beta_pressure_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_beta_pressure, only: Abstract_Beta_Pressure, Concrete_Beta_Pressure, Null_Beta_Pressure

implicit none

private
public :: create, destroy

contains

    subroutine create(beta_pressure, box_size_can_change, generating_data, prefix)
        class(Abstract_Beta_Pressure), allocatable, intent(out) :: beta_pressure
        logical, intent(in) :: box_size_can_change
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: density, excess

        if (box_size_can_change) then
            data_field = prefix//"Beta Pressure.density"
            call generating_data%get(data_field, density, data_found)
            call check_data_found(data_field, data_found)
            data_field = prefix//"Beta Pressure.excess"
            call generating_data%get(data_field, excess, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Beta_Pressure :: beta_pressure)
        else
            allocate(Null_Beta_Pressure :: beta_pressure)
        end if
        call beta_pressure%set(density, excess)
    end subroutine create

    subroutine destroy(beta_pressure)
        class(Abstract_Beta_Pressure), allocatable, intent(inout) :: beta_pressure

        if (allocated(beta_pressure)) deallocate(beta_pressure)
    end subroutine destroy

end module procedures_beta_pressure_factory
