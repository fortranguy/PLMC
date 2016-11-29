module procedures_beta_pressures_excess_factory

use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_beta_pressure_excess, only: Abstract_Beta_Pressure_Excess, XYZ_Beta_Pressure_Excess, &
    XY_Beta_Pressure_Excess, Null_Beta_Pressure_Excess
use procedures_environment_inquirers, only: periodicity_is_xyz, periodicity_is_xy

implicit none

private
public :: create, destroy

contains

    subroutine create(beta_pressures_excess, periodic_boxes, accessible_domains, count_contacts)
        class(Abstract_Beta_Pressure_Excess), allocatable, intent(out) :: beta_pressures_excess(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        class(Abstract_Parallelepiped_Domain), intent(in) :: accessible_domains(:)
        logical, intent(in) :: count_contacts

        integer :: i_box

        if (count_contacts) then
            if (all(periodicity_is_xyz(periodic_boxes))) then
                allocate(XYZ_Beta_Pressure_Excess :: beta_pressures_excess(size(periodic_boxes)))
            else if (all(periodicity_is_xy(periodic_boxes))) then
                allocate(XY_Beta_Pressure_Excess :: beta_pressures_excess(size(periodic_boxes)))
            else
                call error_exit("procedures_beta_pressures_excess_factory: create: "//&
                        "box periodicity is unknown.")
            end if
        else
            allocate(Null_Beta_Pressure_Excess :: beta_pressures_excess(size(periodic_boxes)))
        end if

        do i_box = 1, size(beta_pressures_excess)
            call beta_pressures_excess(i_box)%construct(accessible_domains(i_box))
        end do
    end subroutine create

    subroutine destroy(beta_pressures_excess)
        class(Abstract_Beta_Pressure_Excess), allocatable, intent(inout) :: beta_pressures_excess(:)

        integer :: i_box

        if (allocated(beta_pressures_excess)) then
            do i_box = size(beta_pressures_excess), 1, -1
                call beta_pressures_excess(i_box)%destroy()
            end do
            deallocate(beta_pressures_excess)
        end if
    end subroutine destroy

end module procedures_beta_pressures_excess_factory
