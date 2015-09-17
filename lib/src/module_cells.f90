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
public Abstract_Cells_construct, Abstract_Cells_destroy

    integer, parameter :: num_local_cells(num_dimensions) = 3

    type, public :: Abstract_Cells
        integer :: nums(num_dimensions)
        real(DP) :: size(num_dimensions)
        integer, dimension(num_dimensions) :: global_lbounds, global_ubounds
        integer, dimension(num_dimensions) :: local_lbounds, local_ubounds
        class(Abstract_Visitable_List), allocatable :: visitables_lists(:, :, :)
        integer, allocatable :: neighbours(:, :, :, :, :, :, :)
    end type Abstract_Cells

    type, extends(Abstract_Cells), public :: XYZ_Cells

    end type XYZ_Cells

contains

    subroutine Abstract_Cells_construct(cells, mold, periodic_box, min_cell_edge)
        type(Abstract_Cells), intent(out) :: cells
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Visitable_List), intent(in) :: mold
        real(DP), intent(in) :: min_cell_edge

        cells%nums = floor(periodic_box%get_size()/min_cell_edge)
        call check_positive("Abstract_Cells", "cells%nums", cells%nums)
        cells%size = periodic_box%get_size() / real(cells%nums, DP)
        call check_positive("Abstract_Cells", "cells%size", cells%size)
        call Abstract_Cells_check_size(cells, periodic_box)

        cells%global_lbounds = -cells%nums/2
        cells%global_ubounds = cells%global_lbounds + cells%nums - 1
        allocate(cells%visitables_lists(cells%global_lbounds(1):cells%global_ubounds(1), &
                                         cells%global_lbounds(2):cells%global_ubounds(2), &
                                         cells%global_lbounds(3):cells%global_ubounds(3)), &
                                         mold=mold)
        cells%local_lbounds = -num_local_cells/2
        cells%local_ubounds = cells%local_lbounds + num_local_cells - 1
        allocate(cells%neighbours(3, cells%local_lbounds(1):cells%local_ubounds(1), &
                                     cells%local_lbounds(2):cells%local_ubounds(2), &
                                     cells%local_lbounds(3):cells%local_ubounds(3), &
                                     cells%global_lbounds(1):cells%global_ubounds(1), &
                                     cells%global_lbounds(2):cells%global_ubounds(2), &
                                     cells%global_lbounds(3):cells%global_ubounds(3)))
    end subroutine Abstract_Cells_construct

    subroutine Abstract_Cells_check_size(cells, periodic_box)
        type(Abstract_Cells), intent(in) :: cells
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        real(DP) :: box_mod_cell(num_dimensions)

        if (any(cells%nums < num_local_cells)) then
            call error_exit("Abstract_Cells: cells%nums is too small.")
        end if
        box_mod_cell = modulo(periodic_box%get_size(), cells%size)
        if (any(box_mod_cell > real_zero .and. abs(box_mod_cell - cells%size) > real_zero)) then
            call error_exit("Abstract_Cells:"//&
                            "cells%size size is not a divisor of periodic_box%get_size()")
        end if
    end subroutine Abstract_Cells_check_size

    subroutine Abstract_Cells_destroy(cells)
        type(Abstract_Cells), intent(inout) :: cells

        if (allocated(cells%neighbours)) deallocate(cells%neighbours)
        if (allocated(cells%visitables_lists)) deallocate(cells%visitables_lists)
    end subroutine Abstract_Cells_destroy

    subroutine XYZ_Cells_set_neighbours(cells)
        type(Abstract_Cells), intent(inout) :: cells

        integer :: global_i1, global_i2, global_i3
        integer :: local_i1, local_i2, local_i3
        integer :: i_cell(num_dimensions)

        cells%neighbours = 0
        do global_i3 = cells%global_lbounds(3), cells%global_ubounds(3)
        do global_i2 = cells%global_lbounds(2), cells%global_ubounds(2)
        do global_i1 = cells%global_lbounds(1), cells%global_ubounds(1)
            do local_i3 = cells%local_lbounds(3), cells%local_ubounds(3)
            do local_i2 = cells%local_lbounds(2), cells%local_ubounds(2)
            do local_i1 = cells%local_lbounds(1), cells%local_ubounds(1)
                i_cell = [global_i1, global_i2, global_i3] + [local_i1, local_i2, local_i3]
                i_cell = periodic_i_cell(i_cell, cells%global_lbounds, cells%global_ubounds)
                cells%neighbours(:, local_i1, local_i2, local_i3, &
                    global_i1, global_i2, global_i3) = i_cell
            end do
            end do
            end do
        end do
        end do
        end do
    end subroutine XYZ_Cells_set_neighbours

    pure function periodic_i_cell(i_cell, lbounds, ubounds)
        integer, intent(in) :: i_cell(:), lbounds(:), ubounds(:)
        integer :: periodic_i_cell(3)

        integer :: nums_cells(3)

        nums_cells = ubounds + 1 - lbounds
        periodic_i_cell = i_cell
        where (periodic_i_cell < lbounds)
            periodic_i_cell = periodic_i_cell + nums_cells
        end where
        where (periodic_i_cell > ubounds)
            periodic_i_cell = periodic_i_cell - nums_cells
        end where
    end function periodic_i_cell

end module module_cells
