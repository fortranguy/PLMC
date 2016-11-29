module procedures_dlc_weight_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use classes_permittivity, only: Abstract_Permittivity
use classes_dlc_weight, only: Abstract_DLC_Weight, Concrete_DLC_Weight, Null_DLC_Weight
use procedures_environment_inquirers, only: periodicity_is_xy

implicit none

private
public :: create, destroy

contains

    subroutine create(weight, periodic_box, reciprocal_lattice, permittivity, dipoles_exist)
        class(Abstract_DLC_Weight), allocatable, intent(out) :: weight
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Permittivity), intent(in) :: permittivity
        logical, intent(in) :: dipoles_exist

        if (periodicity_is_xy(periodic_box) .and. dipoles_exist) then
            allocate(Concrete_DLC_Weight :: weight)
        else
            allocate(Null_DLC_Weight :: weight)
        end if
        call weight%construct(periodic_box, reciprocal_lattice, permittivity)
    end subroutine create

    subroutine destroy(weight)
        class(Abstract_DLC_Weight), allocatable, intent(inout) :: weight

        if (allocated(weight)) then
            call weight%destroy()
            deallocate(weight)
        end if
    end subroutine destroy

end module procedures_dlc_weight_factory
