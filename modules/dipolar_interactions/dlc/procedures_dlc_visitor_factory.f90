module procedures_dlc_visitor_factory

use procedures_errors, only: error_exit
use types_environment_wrapper, only: Environment_Wrapper
use classes_dlc_weight, only: Abstract_DLC_Weight
use classes_dlc_structures, only: Abstract_DLC_Structures, Concrete_DLC_Structures, &
    Null_DLC_Structures
use classes_dlc_visitor, only: Abstract_DLC_Visitor, Concrete_DLC_Visitor, Null_DLC_Visitor

implicit none

private
public :: create, destroy

contains

    subroutine create(visitor, environment, weight, structures)
        class(Abstract_DLC_Visitor), allocatable, intent(out) :: visitor
        type(Environment_Wrapper), intent(in) :: environment
        class(Abstract_DLC_Weight), intent(in) :: weight
        class(Abstract_DLC_Structures), intent(in) :: structures

        select type (structures)
            type is (Concrete_DLC_Structures)
                allocate(Concrete_DLC_Visitor :: visitor)
            type is (Null_DLC_Structures)
                allocate(Null_DLC_Visitor :: visitor)
            class default
                call error_exit("procedures_dlc_visitor_factory: create: structures type unknown.")
        end select
        call visitor%construct(environment%periodic_box, environment%reciprocal_lattice, weight, &
            structures)
    end subroutine create

    subroutine destroy(visitor)
        class(Abstract_DLC_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy

end module procedures_dlc_visitor_factory
