module module_cells

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use data_geometry, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive
use class_periodic_box, only: Abstract_Periodic_Box
use class_positions, only: Abstract_Positions
use class_pair_potential, only: Abstract_Pair_Potential
use class_visitable_list, only: Abstract_Visitable_List

implicit none

private
public Concrete_Cells_construct, Concrete_Cells_destroy

    integer, parameter :: num_local_cells(num_dimensions) = 3

    type, public :: Concrete_Cells
        integer :: nums(num_dimensions)
        real(DP) :: size(num_dimensions)
        integer, dimension(num_dimensions) :: global_lbounds, global_ubounds
        integer, dimension(num_dimensions) :: local_lbounds, local_ubounds
        class(Abstract_Visitable_List), allocatable :: visitables_lists(:, :, :)
        integer, allocatable :: neighbours(:, :, :, :, :, :, :)
        procedure(Concrete_Cells_set_neighbours), pointer, private :: set_neighbours
    contains
        procedure, private :: check_size => Concrete_Cells_check_size
        procedure, private :: construct_visitables_lists => &
            Concrete_Cells_construct_visitables_lists
    end type Concrete_Cells

    abstract interface

        subroutine Concrete_Cells_set_neighbours(this)
        import :: Concrete_Cells
            class(Concrete_Cells), intent(inout) :: this
        end subroutine Concrete_Cells_set_neighbours

    end interface

contains

    subroutine Concrete_Cells_construct(this, mold, periodic_box, periodicity, positions, &
            pair_potential)
        class(Concrete_Cells), intent(out) :: this
        class(Abstract_Visitable_List), intent(in) :: mold
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        character(len=*), intent(in) :: periodicity
        class(Abstract_Positions), intent(in) :: positions
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        this%nums = floor(periodic_box%get_size()/pair_potential%get_max_distance())
        call check_positive("Concrete_Cells", "this%nums", this%nums)
        this%size = periodic_box%get_size() / real(this%nums, DP)
        call check_positive("Concrete_Cells", "this%size", this%size)
        call this%check_size(periodic_box)

        this%global_lbounds = -this%nums/2
        this%global_ubounds = this%global_lbounds + this%nums - 1
        allocate(this%visitables_lists(this%global_lbounds(1):this%global_ubounds(1), &
                                         this%global_lbounds(2):this%global_ubounds(2), &
                                         this%global_lbounds(3):this%global_ubounds(3)), &
                                         mold=mold)
        call this%construct_visitables_lists(periodic_box, positions)

        this%local_lbounds = -num_local_cells/2
        this%local_ubounds = this%local_lbounds + num_local_cells - 1
        allocate(this%neighbours(3, this%local_lbounds(1):this%local_ubounds(1), &
                                     this%local_lbounds(2):this%local_ubounds(2), &
                                     this%local_lbounds(3):this%local_ubounds(3), &
                                     this%global_lbounds(1):this%global_ubounds(1), &
                                     this%global_lbounds(2):this%global_ubounds(2), &
                                     this%global_lbounds(3):this%global_ubounds(3)))

        select case (periodicity)
            case ("3D")
                this%set_neighbours => Concrete_Cells_set_neighbours_PBC_3D
            case default
                call error_exit(periodicity//" unknown.")
        end select
        call this%set_neighbours()
    end subroutine Concrete_Cells_construct

    subroutine Concrete_Cells_construct_visitables_lists(this, periodic_box, positions)
        class(Concrete_Cells), intent(inout) :: this
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Positions), intent(in) :: positions

        integer :: global_i1, global_i2, global_i3

        do global_i3 = this%global_lbounds(3), this%global_ubounds(3)
        do global_i2 = this%global_lbounds(2), this%global_ubounds(2)
        do global_i1 = this%global_lbounds(1), this%global_ubounds(1)
            call this%visitables_lists(global_i1, global_i2, global_i3)%construct(periodic_box, &
                positions)
        end do
        end do
        end do
    end subroutine Concrete_Cells_construct_visitables_lists

    subroutine Concrete_Cells_check_size(this, periodic_box)
        class(Concrete_Cells), intent(in) :: this
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        real(DP) :: box_mod_cell(num_dimensions)

        if (any(this%nums < num_local_cells)) then
            call error_exit("Concrete_Cells: this%nums is too small.")
        end if
        box_mod_cell = modulo(periodic_box%get_size(), this%size)
        if (any(box_mod_cell > real_zero .and. abs(box_mod_cell - this%size) > real_zero)) then
            call error_exit("Concrete_Cells:"//&
                            "this%size size is not a divisor of periodic_box%get_size()")
        end if
    end subroutine Concrete_Cells_check_size

    subroutine Concrete_Cells_destroy(this)
        class(Concrete_Cells), intent(inout) :: this

        integer :: global_i1, global_i2, global_i3

        do global_i3 = this%global_ubounds(3), this%global_lbounds(3), -1
        do global_i2 = this%global_ubounds(2), this%global_lbounds(2), -1
        do global_i1 = this%global_ubounds(1), this%global_lbounds(1), -1
            call this%visitables_lists(global_i1, global_i2, global_i3)%destroy()
        end do
        end do
        end do

        if (allocated(this%neighbours)) deallocate(this%neighbours)
        if (allocated(this%visitables_lists)) deallocate(this%visitables_lists)
    end subroutine Concrete_Cells_destroy

    subroutine Concrete_Cells_set_neighbours_PBC_3D(this)
        class(Concrete_Cells), intent(inout) :: this

        integer :: global_i1, global_i2, global_i3
        integer :: local_i1, local_i2, local_i3
        integer :: i_cell(num_dimensions)

        this%neighbours = 0
        do global_i3 = this%global_lbounds(3), this%global_ubounds(3)
        do global_i2 = this%global_lbounds(2), this%global_ubounds(2)
        do global_i1 = this%global_lbounds(1), this%global_ubounds(1)
            do local_i3 = this%local_lbounds(3), this%local_ubounds(3)
            do local_i2 = this%local_lbounds(2), this%local_ubounds(2)
            do local_i1 = this%local_lbounds(1), this%local_ubounds(1)
                i_cell = [global_i1, global_i2, global_i3] + [local_i1, local_i2, local_i3]
                i_cell = pbc_3d(i_cell, this%global_lbounds, this%global_ubounds)
                this%neighbours(:, local_i1, local_i2, local_i3, &
                    global_i1, global_i2, global_i3) = i_cell
            end do
            end do
            end do
        end do
        end do
        end do
    end subroutine Concrete_Cells_set_neighbours_PBC_3D

    pure function pbc_3d(i_cell, lbounds, ubounds)
        integer, intent(in) :: i_cell(:), lbounds(:), ubounds(:)
        integer :: pbc_3d(3)

        integer :: nums_this(3)

        nums_this = ubounds + 1 - lbounds
        pbc_3d = i_cell
        where (pbc_3d < lbounds)
            pbc_3d = pbc_3d + nums_this
        end where
        where (pbc_3d > ubounds)
            pbc_3d = pbc_3d - nums_this
        end where
    end function pbc_3d

end module module_cells
