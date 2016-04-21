module procedures_temperature_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_temperature, only: Abstract_Temperature, Concrete_Temperature

implicit none

private
public :: create, destroy

contains

    subroutine create(temperature, input_data, prefix)
        class(Abstract_Temperature), allocatable, intent(out) :: temperature
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        real(DP) :: temperature_value
        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = prefix//"Thermostat.temperature"
        call input_data%get(data_field, temperature_value, data_found)
        call check_data_found(data_field, data_found)
        allocate(Concrete_Temperature :: temperature)
        call temperature%set(temperature_value)
    end subroutine create

    subroutine destroy(temperature)
        class(Abstract_Temperature), allocatable, intent(inout) :: temperature

        if (allocated(temperature)) deallocate(temperature)
    end subroutine destroy

end module procedures_temperature_factory
