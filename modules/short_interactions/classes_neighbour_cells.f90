module classes_neighbour_cells

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, real_zero
use data_cells, only: nums_local_cells
use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_pair_potential, only: Abstract_Pair_Potential
use procedures_neighbour_cells_micro, only: pbc_3d_index

implicit none

private

    type, abstract, public :: Abstract_Neighbour_Cells
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Pair_Potential), pointer :: pair_potential => null()
        integer :: nums(num_dimensions) = 0
        real(DP) :: size(num_dimensions) = 0._DP
        integer, dimension(num_dimensions) :: global_lbounds = 0, global_ubounds = 0
        logical :: skip_bottom_layer(-nums_local_cells(3)/2:nums_local_cells(3)/2) = .false.
        logical :: skip_top_layer(-nums_local_cells(3)/2:nums_local_cells(3)/2) = .false.
        integer, allocatable :: neighbours(:, :, :, :, :, :, :)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset => Abstract_reset
        procedure :: get_global_lbounds => Abstract_get_global_lbounds
        procedure :: get_global_ubounds => Abstract_get_global_ubounds
        procedure :: get => Abstract_get
        procedure :: skip => Abstract_skip
        procedure :: index => Abstract_index
        procedure(Abstract_set_skip_layers), private, deferred :: set_skip_layers
        procedure(Abstract_check_nums), private, deferred :: check_nums
        procedure, private :: check_size => Abstract_check_size
        procedure, private :: set_neighbours => Abstract_set_neighbours
    end type Abstract_Neighbour_Cells

    abstract interface

        pure subroutine Abstract_set_skip_layers(this)
        import :: Abstract_Neighbour_Cells
            class(Abstract_Neighbour_Cells), intent(inout) :: this
        end subroutine Abstract_set_skip_layers

        subroutine Abstract_check_nums(this)
        import :: Abstract_Neighbour_Cells
            class(Abstract_Neighbour_Cells), intent(in) :: this
        end subroutine Abstract_check_nums

    end interface

    type, extends(Abstract_Neighbour_Cells), public :: XYZ_PBC_Neighbour_Cells
    contains
        procedure, private :: check_nums => XYZ_check_nums
        procedure, private :: set_skip_layers => XYZ_set_skip_layers
    end type XYZ_PBC_Neighbour_Cells

    type, extends(Abstract_Neighbour_Cells), public :: XY_PBC_Neighbour_Cells
    contains
        procedure, private :: check_nums => XY_check_nums
        procedure, private :: set_skip_layers => XY_set_skip_layers
    end type XY_PBC_Neighbour_Cells

contains

!implementation Abstract_Neighbour_Cells

    subroutine Abstract_construct(this, periodic_box, pair_potential)
        class(Abstract_Neighbour_Cells), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential

        this%periodic_box => periodic_box
        this%pair_potential => pair_potential
        call this%set_skip_layers()
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Neighbour_Cells), intent(inout) :: this

        if (allocated(this%neighbours)) deallocate(this%neighbours)
        this%pair_potential => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_reset(this)
        class(Abstract_Neighbour_Cells), intent(inout) :: this

        this%nums = floor(this%periodic_box%get_size()/this%pair_potential%get_max_distance())
        call this%check_nums()
        this%size = this%periodic_box%get_size() / real(this%nums, DP)
        call this%check_size()

        this%global_lbounds = -this%nums/2
        this%global_ubounds = this%global_lbounds + this%nums - 1

        if (allocated(this%neighbours)) deallocate(this%neighbours)
        allocate(this%neighbours(3, -nums_local_cells(1)/2:nums_local_cells(1)/2, &
                                    -nums_local_cells(2)/2:nums_local_cells(2)/2, &
                                    -nums_local_cells(3)/2:nums_local_cells(3)/2, &
                                    this%global_lbounds(1):this%global_ubounds(1), &
                                    this%global_lbounds(2):this%global_ubounds(2), &
                                    this%global_lbounds(3):this%global_ubounds(3)))
        call this%set_neighbours()
    end subroutine Abstract_reset

    subroutine Abstract_check_size(this)
        class(Abstract_Neighbour_Cells), intent(in) :: this

        real(DP) :: box_modulo_cell(num_dimensions)

        box_modulo_cell = modulo(this%periodic_box%get_size(), this%size)
        if (any(box_modulo_cell > real_zero .and. abs(box_modulo_cell-this%size) > real_zero)) then
            call error_exit("Abstract_Neighbour_Cells: check_size: "//&
                            "this%size is not a divisor of periodic_box%get_size()")
        end if
    end subroutine Abstract_check_size

    subroutine Abstract_set_neighbours(this)
        class(Abstract_Neighbour_Cells), intent(inout) :: this

        integer :: global_i1, global_i2, global_i3
        logical :: at_bottom_layer, at_top_layer
        integer :: local_i1, local_i2, local_i3
        integer :: ijk_cell(num_dimensions)

        this%neighbours = 0
        do global_i3 = this%global_lbounds(3), this%global_ubounds(3)
            at_bottom_layer = (global_i3 == this%global_lbounds(3))
            at_top_layer = (global_i3 == this%global_ubounds(3))
        do global_i2 = this%global_lbounds(2), this%global_ubounds(2)
        do global_i1 = this%global_lbounds(1), this%global_ubounds(1)
            do local_i3 = -nums_local_cells(3)/2, nums_local_cells(3)/2
                if (this%skip(at_bottom_layer, at_top_layer, local_i3)) cycle
            do local_i2 = -nums_local_cells(2)/2, nums_local_cells(2)/2
            do local_i1 = -nums_local_cells(1)/2, nums_local_cells(1)/2
                ijk_cell = [global_i1, global_i2, global_i3] + [local_i1, local_i2, local_i3]
                ijk_cell = pbc_3d_index(ijk_cell, this%nums)
                this%neighbours(:, local_i1, local_i2, local_i3, global_i1, global_i2, global_i3) =&
                    ijk_cell
            end do
            end do
            end do
        end do
        end do
        end do
    end subroutine Abstract_set_neighbours

    pure function Abstract_get_global_lbounds(this) result(lbounds)
        class(Abstract_Neighbour_Cells), intent(in) :: this
        integer :: lbounds(num_dimensions)

        lbounds = this%global_lbounds
    end function Abstract_get_global_lbounds

    pure function Abstract_get_global_ubounds(this) result(ubounds)
        class(Abstract_Neighbour_Cells), intent(in) :: this
        integer :: ubounds(num_dimensions)

        ubounds = this%global_ubounds
    end function Abstract_get_global_ubounds

    pure function Abstract_get(this, local_i1, local_i2, local_i3, global_i1, global_i2, global_i3)&
        result(neighbour)
        class(Abstract_Neighbour_Cells), intent(in) :: this
        integer, intent(in) :: local_i1, local_i2, local_i3
        integer, intent(in) :: global_i1, global_i2, global_i3
        integer :: neighbour(num_dimensions)

        neighbour = this%neighbours(:, local_i1, local_i2, local_i3, global_i1, global_i2, &
            global_i3)
    end function Abstract_get

    pure logical function Abstract_skip(this, at_bottom_layer, at_top_layer, local_i3) result(skip)
        class(Abstract_Neighbour_Cells), intent(in) :: this
        logical, intent(in) :: at_bottom_layer, at_top_layer
        integer, intent(in) :: local_i3

        skip = (at_bottom_layer .and. this%skip_bottom_layer(local_i3)) .or. &
            (at_top_layer .and. this%skip_top_layer(local_i3))
    end function Abstract_skip

    pure function Abstract_index(this, position) result(index)
        class(Abstract_Neighbour_Cells), intent(in) :: this
        real(DP), intent(in) :: position(:)
        integer :: index(num_dimensions)

        where (mod(this%nums, 2) == 0)
            index = floor(position/this%size)
        elsewhere
            index = nint(position/this%size)
        end where
    end function Abstract_index

!end implementation Abstract_Neighbour_Cells


!implementation XYZ_PBC_Neighbour_Cells

    pure subroutine XYZ_set_skip_layers(this)
        class(XYZ_PBC_Neighbour_Cells), intent(inout) :: this

        this%skip_top_layer = .false.
        this%skip_bottom_layer = .false.
    end subroutine XYZ_set_skip_layers

    subroutine XYZ_check_nums(this)
        class(XYZ_PBC_Neighbour_Cells), intent(in) :: this

        if (any(this%nums < nums_local_cells)) then
            call error_exit("XYZ_PBC_Neighbour_Cells: this%nums is too small.")
        end if
    end subroutine XYZ_check_nums

!end implementation XYZ_PBC_Neighbour_Cells

!implementation XY_PBC_Neighbour_Cells

    pure subroutine XY_set_skip_layers(this)
        class(XY_PBC_Neighbour_Cells), intent(inout) :: this

        integer :: local_i3

        do local_i3 = -nums_local_cells(3)/2, nums_local_cells(3)/2
            if (local_i3 == 1) then
                this%skip_top_layer(local_i3) = .true.
            else
                this%skip_top_layer(local_i3) = .false.
            end if
            if (local_i3 == -1) then
                this%skip_bottom_layer(local_i3) = .true.
            else
                this%skip_bottom_layer(local_i3) = .false.
            end if
        end do
    end subroutine XY_set_skip_layers

    subroutine XY_check_nums(this)
        class(XY_PBC_Neighbour_Cells), intent(in) :: this

        if (any(this%nums(1:2) < nums_local_cells(1:2))) then
            call error_exit("XY_PBC_Neighbour_Cells: this%nums is too small.")
        end if
    end subroutine XY_check_nums

!end implementation XY_PBC_Neighbour_Cells

end module classes_neighbour_cells
