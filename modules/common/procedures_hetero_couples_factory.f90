module procedures_hetero_couples_factory

use classes_hetero_couples, only: Abstract_Hetero_Couples, Half_Hetero_Couples, &
    Full_Hetero_Couples, Null_Hetero_Couples

implicit none

private
public :: create_half, create_full, destroy

contains

    subroutine create_half(couples, num_partners)
        class(Abstract_Hetero_Couples), allocatable, intent(out) :: couples
        integer, intent(in) :: num_partners

        if (num_partners > 1) then
            allocate(Half_Hetero_Couples :: couples)
        else
            allocate(Null_Hetero_Couples :: couples)
        end if
        call couples%construct(num_partners)
    end subroutine create_half

    subroutine create_full(couples, num_partners)
        class(Abstract_Hetero_Couples), allocatable, intent(out) :: couples
        integer, intent(in) :: num_partners

        if (num_partners > 1) then
            allocate(Full_Hetero_Couples :: couples)
        else
            allocate(Null_Hetero_Couples :: couples)
        end if
        call couples%construct(num_partners)
    end subroutine create_full

    subroutine destroy(couples)
        class(Abstract_Hetero_Couples), allocatable, intent(inout) :: couples

        if (allocated(couples)) then
            call couples%destroy()
            deallocate(couples)
        end if
    end subroutine destroy

end module procedures_hetero_couples_factory
