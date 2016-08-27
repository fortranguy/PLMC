module procedures_beta_pressure_excess_factory

use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_beta_pressure_excess, only: Abstract_Beta_Pressure_Excess, XYZ_Beta_Pressure_Excess, &
    XYZ_Beta_Pressure_Excess, Null_Beta_Pressure_Excess
use procedures_environment_inquirers, only: periodicity_is_xyz, periodicity_is_xy

implicit none

private
public :: create, destroy

contains

    subroutine create(beta_pressure_excess, periodic_box, accessible_domain, count_contacts)
        class(Abstract_Beta_Pressure_Excess), allocatable, intent(out) :: beta_pressure_excess
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Parallelepiped_Domain), intent(in) :: accessible_domain
        logical, intent(in) :: count_contacts

        if (count_contacts) then
            if (periodicity_is_xyz(periodic_box)) then
                allocate(XYZ_Beta_Pressure_Excess :: beta_pressure_excess)
            else if (periodicity_is_xy(periodic_box)) then
                allocate(XYZ_Beta_Pressure_Excess :: beta_pressure_excess)
            else
                call error_exit("procedures_beta_pressure_excess_factory: create: "//&
                        "box periodicity is unknown.")
            end if
        else
            allocate(Null_Beta_Pressure_Excess :: beta_pressure_excess)
        end if
        call beta_pressure_excess%construct(accessible_domain)
    end subroutine create

    subroutine destroy(beta_pressure_excess)
        class(Abstract_Beta_Pressure_Excess), allocatable, intent(inout) :: beta_pressure_excess

        if (allocated(beta_pressure_excess)) then
            call beta_pressure_excess%destroy()
            deallocate(beta_pressure_excess)
        end if
    end subroutine destroy

end module procedures_beta_pressure_excess_factory
