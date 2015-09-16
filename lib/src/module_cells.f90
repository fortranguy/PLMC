module module_cells

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use data_geometry, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive
use class_periodic_box, only: Abstract_Periodic_Box
use class_visitable_list, only: Abstract_Visitable_List

implicit none

private
public Concrete_Cells_construct, Concrete_Cells_destroy

    integer, parameter :: nums_neighbour_cells(num_dimensions) = 3

    type, public :: Concrete_Cells
        integer :: nums(num_dimensions)
        real(DP) :: size(num_dimensions)
        class(Abstract_Visitable_List), allocatable :: visitables_lists(:, :, :)
    end type Concrete_Cells

contains

    subroutine Concrete_Cells_construct(cells, mold, periodic_box, min_cell_edge)
        type(Concrete_Cells), intent(out) :: cells
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Visitable_List), intent(in) :: mold
        real(DP), intent(in) :: min_cell_edge

        integer, dimension(num_dimensions) :: mi_bounds

        cells%nums = floor(periodic_box%get_size()/min_cell_edge)
        call check_positive("Concrete_Cells", "cells%nums", cells%nums)
        cells%size = periodic_box%get_size() / real(cells%nums, DP)
        call check_positive("Concrete_Cells", "cells%size", cells%size)
        call Concrete_Cells_check_size(cells, periodic_box)
        mi_bounds = int(real(cells%nums, DP)/2._DP)
        allocate(cells%visitables_lists(-mi_bounds(1):mi_bounds(1), -mi_bounds(2):mi_bounds(2), &
                                        -mi_bounds(3):mi_bounds(3)), mold=mold)
    end subroutine Concrete_Cells_construct

    subroutine Concrete_Cells_check_size(cells, periodic_box)
        type(Concrete_Cells), intent(in) :: cells
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        real(DP) :: box_mod_cell(num_dimensions)

        if (any(cells%nums < nums_neighbour_cells)) then
            call error_exit("Concrete_Cells: cells%nums is too small.")
        end if
        box_mod_cell = modulo(periodic_box%get_size(), cells%size)
        if (any(box_mod_cell > real_zero .and. abs(box_mod_cell - cells%size) > real_zero)) then
            call error_exit("Concrete_Cells:"//&
                            "cells%size size is not a divisor of periodic_box%get_size()")
        end if
    end subroutine Concrete_Cells_check_size

    subroutine Concrete_Cells_destroy(cells)
        type(Concrete_Cells), intent(inout) :: cells

        if (allocated(cells%visitables_lists)) deallocate(cells%visitables_lists)
    end subroutine Concrete_Cells_destroy

end module module_cells
