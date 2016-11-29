module procedures_hard_contact_factory

use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_environment_inquirers, only: periodicity_is_xyz, periodicity_is_xy
use classes_dirac_distribution_plus, only: Abstract_Dirac_Distribution_Plus
use classes_hard_contact, only: Abstract_Hard_Contact, XYZ_Hard_Contact, XY_Hard_Contact, &
    Null_Hard_Contact

implicit none

private
public :: create, destroy

contains

    subroutine create(hard_contact, periodic_boxes, dirac_plus, count_contacts)
        class(Abstract_Hard_Contact), allocatable, intent(out) :: hard_contact
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        class(Abstract_Dirac_Distribution_Plus), intent(in) :: dirac_plus
        logical, intent(in) :: count_contacts

        if (count_contacts) then
            if (all(periodicity_is_xyz(periodic_boxes))) then
                allocate(XYZ_Hard_Contact :: hard_contact)
            else if (all(periodicity_is_xy(periodic_boxes))) then
                allocate(XY_Hard_Contact :: hard_contact)
            else
                call error_exit("procedures_hard_contact_factory: create: "//&
                        "box periodicity is unknown.")
            end if
        else
            allocate(Null_Hard_Contact :: hard_contact)
        end if
        call hard_contact%construct(dirac_plus)
    end subroutine create

    subroutine destroy(hard_contact)
        class(Abstract_Hard_Contact), allocatable, intent(inout) :: hard_contact

        if (allocated(hard_contact)) then
            call hard_contact%destroy()
            deallocate(hard_contact)
        end if
    end subroutine destroy

end module procedures_hard_contact_factory
