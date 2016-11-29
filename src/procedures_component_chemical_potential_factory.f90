module procedures_component_chemical_potential_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_component_chemical_potential, only : Abstract_Component_Chemical_Potential, &
    Concrete_Component_Chemical_Potential, Null_Component_Chemical_Potential

implicit none

private
public :: create, destroy

contains

    subroutine create(chemical_potential, can_exchange, generating_data, prefix)
        class(Abstract_Component_Chemical_Potential), allocatable, intent(out) :: &
            chemical_potential
        logical, intent(in) :: can_exchange
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: density, inv_pow_activity

        if (can_exchange) then
            data_field = prefix//"density"
            call generating_data%get(data_field, density, data_found)
            call check_data_found(data_field, data_found)
            data_field = prefix//"inverse powered activity"
            call generating_data%get(data_field, inv_pow_activity, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Component_Chemical_Potential :: chemical_potential)
        else
            allocate(Null_Component_Chemical_Potential :: chemical_potential)
        end if
        call chemical_potential%set(density, inv_pow_activity)
    end subroutine create

    subroutine destroy(chemical_potential)
        class(Abstract_Component_Chemical_Potential), allocatable, intent(inout) :: &
            chemical_potential

        if (allocated(chemical_potential)) deallocate(chemical_potential)
    end subroutine destroy

end module procedures_component_chemical_potential_factory
