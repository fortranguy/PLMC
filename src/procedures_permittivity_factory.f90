module procedures_permittivity_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_permittivity, only: Abstract_Permittivity, Concrete_Permittivity, Null_Permittivity
use procedures_environment_inquirers, only: use_permittivity

implicit none

private
public :: create, destroy

contains

    subroutine create(permittivity, generating_data, prefix)
        class(Abstract_Permittivity), allocatable, intent(out) :: permittivity
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: permittivity_value

        if (use_permittivity(generating_data, prefix)) then
            data_field = prefix//"Permittivity.value"
            call generating_data%get(data_field, permittivity_value, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Permittivity :: permittivity)
        else
            allocate(Null_Permittivity :: permittivity)
        end if
        call permittivity%set(permittivity_value)
    end subroutine create

    subroutine destroy(permittivity)
        class(Abstract_Permittivity), allocatable, intent(inout) :: permittivity

        if (allocated(permittivity)) then
            deallocate(permittivity)
        end if
    end subroutine destroy

end module procedures_permittivity_factory
