module class_component_positions

use, intrinsic :: iso_fortran_env, only: DP => REAL64, error_unit
use data_constants, only: num_dimensions
use procedures_checks, only: check_in_range, check_3d_array
use class_periodic_box, only: Abstract_Periodic_Box
use class_component_number, only: Abstract_Component_Number
use procedures_coordinates_micro, only: increase_coordinates_size
use class_component_coordinates, only: Abstract_Component_Coordinates

implicit none

private

    type, extends(Abstract_Component_Coordinates), abstract, public :: Abstract_Component_Positions
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box
        class(Abstract_Component_Number), pointer :: number
        real(DP), allocatable :: positions(:, :)
    contains
        procedure :: construct => Abstract_Component_Positions_construct
        procedure :: destroy => Abstract_Component_Positions_destroy
        procedure :: set => Abstract_Component_Positions_set
        procedure :: get_num => Abstract_Component_Positions_get_num
        procedure :: get => Abstract_Component_Positions_get
        procedure :: add => Abstract_Component_Positions_add
        procedure :: remove => Abstract_Component_Positions_remove
    end type Abstract_Component_Positions

    type, extends(Abstract_Component_Positions), public :: Concrete_Component_Positions

    end type Concrete_Component_Positions

    type, extends(Abstract_Component_Positions), public :: Null_Component_Positions
    contains
        procedure :: construct => Null_Component_Positions_construct
        procedure :: destroy => Null_Component_Positions_destroy
        procedure :: set => Null_Component_Positions_set
        procedure :: get_num => Null_Component_Positions_get_num
        procedure :: get => Null_Component_Positions_get
        procedure :: add => Null_Component_Positions_add
        procedure :: remove => Null_Component_Positions_remove
    end type Null_Component_Positions

contains

!implementation Abstract_Component_Positions

    subroutine Abstract_Component_Positions_construct(this, periodic_box, number)
        class(Abstract_Component_Positions), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Number), target, intent(in) :: number

        this%periodic_box => periodic_box
        this%number => number
        if (this%number%get() == 0) then
            allocate(this%positions(num_dimensions, 1))
        else
            allocate(this%positions(num_dimensions, this%number%get()))
        end if
        this%positions = 0._DP
    end subroutine Abstract_Component_Positions_construct

    subroutine Abstract_Component_Positions_destroy(this)
        class(Abstract_Component_Positions), intent(inout) :: this

        if (allocated(this%positions)) deallocate(this%positions)
        this%number => null()
        this%periodic_box => null()
    end subroutine Abstract_Component_Positions_destroy

    subroutine Abstract_Component_Positions_set(this, i_particle, vector)
        class(Abstract_Component_Positions), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: vector(:)

        call check_in_range("Abstract_Component_Positions_set", this%number%get(), &
            "i_particle", i_particle)
        call check_3d_array("Abstract_Component_Positions", "vector", vector)
        this%positions(:, i_particle) = this%periodic_box%folded(vector)
    end subroutine Abstract_Component_Positions_set

    pure function Abstract_Component_Positions_get_num(this) result(num_positions)
        class(Abstract_Component_Positions), intent(in) :: this
        integer :: num_positions

        num_positions = this%number%get()
    end function Abstract_Component_Positions_get_num

    pure function Abstract_Component_Positions_get(this, i_particle) result(position)
        class(Abstract_Component_Positions), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: position(num_dimensions)

        position = this%positions(:, i_particle)
    end function Abstract_Component_Positions_get

    subroutine Abstract_Component_Positions_add(this, vector)
        class(Abstract_Component_Positions), intent(inout) :: this
        real(DP), intent(in) :: vector(:)

        if (size(this%positions, 2) < this%number%get()) then
            call increase_coordinates_size(this%positions)
        end if
        call this%set(this%number%get(), vector)
    end subroutine Abstract_Component_Positions_add

    subroutine Abstract_Component_Positions_remove(this, i_particle)
        class(Abstract_Component_Positions), intent(inout) :: this
        integer, intent(in) :: i_particle

        call check_in_range("Abstract_Component_Positions_remove", this%number%get(), &
            "i_particle", i_particle)
        if (i_particle < this%number%get()) then
            call this%set(i_particle, this%get(this%number%get()))
        end if
    end subroutine Abstract_Component_Positions_remove

!end implementation Abstract_Component_Positions

!implementation Null_Component_Positions

    subroutine Null_Component_Positions_construct(this, periodic_box, number)
        class(Null_Component_Positions), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Number), target, intent(in) :: number
    end subroutine Null_Component_Positions_construct

    subroutine Null_Component_Positions_destroy(this)
        class(Null_Component_Positions), intent(inout) :: this
    end subroutine Null_Component_Positions_destroy

    subroutine Null_Component_Positions_set(this, i_particle, vector)
        class(Null_Component_Positions), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: vector(:)
    end subroutine Null_Component_Positions_set

    pure function Null_Component_Positions_get_num(this) result(num_positions)
        class(Null_Component_Positions), intent(in) :: this
        integer :: num_positions
        num_positions = 0
    end function Null_Component_Positions_get_num

    pure function Null_Component_Positions_get(this, i_particle) result(position)
        class(Null_Component_Positions), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: position(num_dimensions)
        position = 0._DP
    end function Null_Component_Positions_get

    subroutine Null_Component_Positions_add(this, vector)
        class(Null_Component_Positions), intent(inout) :: this
        real(DP), intent(in) :: vector(:)
    end subroutine Null_Component_Positions_add

    subroutine Null_Component_Positions_remove(this, i_particle)
        class(Null_Component_Positions), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_Component_Positions_remove

!end implementation Null_Component_Positions

end module class_component_positions
