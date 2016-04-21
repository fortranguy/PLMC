module procedures_des_real_visitor_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_des_real_visitor, only: Abstract_DES_Real_Visitor, Concrete_DES_Real_Visitor, &
    Null_DES_Real_Visitor

implicit none

private
public :: create, destroy

contains

    subroutine create(visitor, periodic_box, dipoles_exist)
        class(Abstract_DES_Real_Visitor), allocatable, intent(out) :: visitor
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: dipoles_exist

        if (dipoles_exist) then
            allocate(Concrete_DES_Real_Visitor :: visitor)
        else
            allocate(Null_DES_Real_Visitor :: visitor)
        end if
        call visitor%construct(periodic_box)
    end subroutine create

    subroutine destroy(visitor)
        class(Abstract_DES_Real_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy

end module procedures_des_real_visitor_factory
