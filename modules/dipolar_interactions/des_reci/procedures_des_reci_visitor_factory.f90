module procedures_des_reci_visitor_factory

use procedures_errors, only: error_exit
use types_environment_wrapper, only: Environment_Wrapper
use classes_des_reci_weight, only: Abstract_DES_Reci_Weight
use classes_des_reci_structure, only: Abstract_DES_Reci_Structure, Concrete_DES_Reci_Structure, &
    Null_DES_Reci_Structure
use classes_des_reci_visitor, only: Abstract_DES_Reci_Visitor, &
    Concrete_DES_Reci_Visitor, Null_DES_Reci_Visitor

implicit none

private
public :: create, destroy

contains

    subroutine create(visitor, environment, weight, structure)
        class(Abstract_DES_Reci_Visitor), allocatable, intent(out) :: visitor
        type(Environment_Wrapper), intent(in) :: environment
        class(Abstract_DES_Reci_Weight), intent(in) :: weight
        class(Abstract_DES_Reci_Structure), intent(in) :: structure

        select type (structure)
            type is (Concrete_DES_Reci_Structure)
                allocate(Concrete_DES_Reci_Visitor :: visitor)
            type is (Null_DES_Reci_Structure)
                allocate(Null_DES_Reci_Visitor :: visitor)
            class default
                call error_exit("procedures_des_reci_visitor_factory: create: structure"//&
                    " type unknown.")
        end select
        call visitor%construct(environment%periodic_box, environment%reciprocal_lattice, &
            weight, structure)
    end subroutine create

    subroutine destroy(visitor)
        class(Abstract_DES_Reci_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy

end module procedures_des_reci_visitor_factory
