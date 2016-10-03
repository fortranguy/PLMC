module procedures_des_reci_weight_factory

use classes_box_size_memento, only: Abstract_Box_Size_Memento
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use classes_permittivity, only: Abstract_Permittivity
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use classes_des_reci_weight, only: Abstract_DES_Reci_Weight, Concrete_DES_Reci_Weight, &
    Null_DES_Reci_Weight

implicit none

private
public :: create, destroy

contains

    subroutine create(weight, box_size_memento, reciprocal_lattice, permittivity, dipoles_exist, &
        alpha)
        class(Abstract_DES_Reci_Weight), allocatable, intent(out) :: weight
        class(Abstract_Box_Size_Memento), intent(in) :: box_size_memento
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Permittivity), intent(in) :: permittivity
        logical, intent(in) :: dipoles_exist
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha

        if (dipoles_exist) then
            allocate(Concrete_DES_Reci_Weight :: weight)
        else
            allocate(Null_DES_Reci_Weight :: weight)
        end if
        call weight%construct(box_size_memento, reciprocal_lattice, permittivity, alpha)
    end subroutine create

    subroutine destroy(weight)
        class(Abstract_DES_Reci_Weight), allocatable, intent(inout) :: weight

        if (allocated(weight)) then
            call weight%destroy()
            deallocate(weight)
        end if
    end subroutine destroy

end module procedures_des_reci_weight_factory
