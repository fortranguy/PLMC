module module_nodes

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions

implicit none

private
public deallocate_list

    type, public :: Concrete_Node
        integer :: i
        real(DP) :: position(num_dimensions)
        type(Concrete_Node), pointer :: next => null()
    end type Concrete_Node

contains

    recursive subroutine deallocate_list(current)
        type(Concrete_Node), pointer, intent(inout) :: current

        if (associated(current%next)) then
            call deallocate_list(current%next)
        end if
        deallocate(current)
    end subroutine deallocate_list

end module module_nodes
