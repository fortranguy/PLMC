module class_neighbour_cells

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, real_zero
use data_cells, only: nums_local_cells
use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Neighbour_Cells
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Pair_Potential), pointer :: pair_potential => null()
        integer :: nums(num_dimensions) = 0
        real(DP) :: size(num_dimensions) = 0._DP
        integer, dimension(num_dimensions) :: lbounds = 0, ubounds = 0
        logical :: skip_bottom_layer(-nums_local_cells(3)/2:nums_local_cells(3)/2) = .false.
        logical :: skip_top_layer(-nums_local_cells(3)/2:nums_local_cells(3)/2) = .false.
        integer, allocatable :: neighbours(:, :, :, :, :, :, :)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: set => Abstract_set
        procedure :: get_lbounds => Abstract_get_lbounds
        procedure :: get_ubounds => Abstract_get_ubounds
        procedure(Abstract_set_skip_layers), private, deferred :: set_skip_layers
        procedure(Abstract_check_nums), private, deferred :: check_nums
        procedure, private :: check_size => Abstract_check_size
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

contains

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

    subroutine Abstract_set(this)
        class(Abstract_Neighbour_Cells), intent(inout) :: this

        this%nums = floor(this%periodic_box%get_size()/this%pair_potential%get_max_distance())
        call this%check_nums()
        this%size = this%periodic_box%get_size() / real(this%nums, DP)
        call this%check_size()

        this%lbounds = -this%nums/2
        this%ubounds = this%lbounds + this%nums - 1

    end subroutine Abstract_set

    subroutine Abstract_check_size(this)
        class(Abstract_Neighbour_Cells), intent(in) :: this

        real(DP) :: box_modulo_cell(num_dimensions)

        box_modulo_cell = modulo(this%periodic_box%get_size(), this%size)
        if (any(box_modulo_cell > real_zero .and. abs(box_modulo_cell-this%size) > real_zero)) then
            call error_exit("Abstract_Neighbour_Cells: check_size: "//&
                            "this%size is not a divisor of periodic_box%get_size()")
        end if
    end subroutine Abstract_check_size

    pure function Abstract_get_lbounds(this) result(lbounds)
        class(Abstract_Neighbour_Cells), intent(in) :: this
        integer :: lbounds(num_dimensions)

        lbounds = this%lbounds
    end function Abstract_get_lbounds

    pure function Abstract_get_ubounds(this) result(ubounds)
        class(Abstract_Neighbour_Cells), intent(in) :: this
        integer :: ubounds(num_dimensions)

        ubounds = this%ubounds
    end function Abstract_get_ubounds

end module class_neighbour_cells
