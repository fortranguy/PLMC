module module_list_node

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions

implicit none

private
public :: deallocate_list, increase_nodes_size

    integer, parameter :: increase_factor = 2

    type :: Concrete_Node
        integer :: i
    end type Concrete_Node

    type, public :: Concrete_Linkable_Node
        integer :: i
        type(Concrete_Linkable_Node), pointer :: next => null()
    end type Concrete_Linkable_Node

contains

    recursive subroutine deallocate_list(current)
        type(Concrete_Linkable_Node), pointer, intent(inout) :: current

        if (associated(current%next)) then
            call deallocate_list(current%next)
        end if
        deallocate(current)
    end subroutine deallocate_list

    subroutine increase_nodes_size(nodes)
        integer, allocatable, intent(inout) :: nodes(:)

        integer, allocatable :: nodes_tmp(:)
        integer :: num_nodes

        num_nodes = size(nodes)
        allocate(nodes_tmp(num_nodes))
        nodes_tmp = nodes
        deallocate(nodes)

        allocate(nodes(increase_factor * num_nodes))
        nodes(1:num_nodes) = nodes_tmp
    end subroutine increase_nodes_size

end module module_list_node
