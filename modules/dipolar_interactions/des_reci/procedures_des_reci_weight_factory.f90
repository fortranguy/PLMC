module procedures_des_reci_weight_factory

use types_environment_wrapper, only: Environment_Wrapper
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use classes_des_reci_weight, only: Abstract_DES_Reci_Weight, Concrete_DES_Reci_Weight, &
    Null_DES_Reci_Weight

implicit none

private
public :: des_reci_weight_create, des_reci_weight_destroy

contains

    subroutine des_reci_weight_create(weight, environment, dipoles_exist, alpha)
        class(Abstract_DES_Reci_Weight), allocatable, intent(out) :: weight
        type(Environment_Wrapper), intent(in) :: environment
        logical, intent(in) :: dipoles_exist
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha

        if (dipoles_exist) then
            allocate(Concrete_DES_Reci_Weight :: weight)
        else
            allocate(Null_DES_Reci_Weight :: weight)
        end if
        call weight%construct(environment%periodic_box, environment%reciprocal_lattice, &
            environment%permittivity, alpha)
    end subroutine des_reci_weight_create

    subroutine des_reci_weight_destroy(weight)
        class(Abstract_DES_Reci_Weight), allocatable, intent(inout) :: weight

        if (allocated(weight)) then
            call weight%destroy()
            deallocate(weight)
        end if
    end subroutine des_reci_weight_destroy

end module procedures_des_reci_weight_factory
