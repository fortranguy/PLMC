module procedures_dlc_weight_factory

use types_environment_wrapper, only: Environment_Wrapper
use classes_dlc_weight, only: Abstract_DLC_Weight, Concrete_DLC_Weight, Null_DLC_Weight
use procedures_environment_inquirers, only: periodicity_is_xy

implicit none

private
public :: create, destroy

contains

    subroutine create(weight, environment, dipoles_exist)
        class(Abstract_DLC_Weight), allocatable, intent(out) :: weight
        type(Environment_Wrapper), intent(in) :: environment
        logical, intent(in) :: dipoles_exist

        if (periodicity_is_xy(environment%periodic_box) .and. dipoles_exist) then
            allocate(Concrete_DLC_Weight :: weight)
        else
            allocate(Null_DLC_Weight :: weight)
        end if
        call weight%construct(environment%periodic_box, environment%reciprocal_lattice, &
            environment%permittivity)
    end subroutine create

    subroutine destroy(weight)
        class(Abstract_DLC_Weight), allocatable, intent(inout) :: weight

        if (allocated(weight)) then
            call weight%destroy()
            deallocate(weight)
        end if
    end subroutine destroy

end module procedures_dlc_weight_factory
