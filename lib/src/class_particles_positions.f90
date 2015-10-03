module class_particles_positions

use, intrinsic :: iso_fortran_env, only: DP => REAL64, error_unit
use data_constants, only: num_dimensions
use procedures_checks, only: check_in_range, check_3d_array
use class_periodic_box, only: Abstract_Periodic_Box
use class_particles_number, only: Abstract_Particles_Number
use procedures_coordinates, only: increase_coordinates_size

implicit none

private

    type, abstract, public :: Abstract_Particles_Positions
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box
        class(Abstract_Particles_Number), pointer :: number
        real(DP), allocatable :: positions(:, :)
    contains
        procedure :: construct => Abstract_Particles_Positions_construct
        procedure :: destroy => Abstract_Particles_Positions_destroy
        procedure :: set => Abstract_Particles_Positions_set
        procedure :: get_num => Abstract_Particles_Positions_get_num
        procedure :: get => Abstract_Particles_Positions_get
        procedure :: add => Abstract_Particles_Positions_add
        procedure :: remove => Abstract_Particles_Positions_remove
    end type Abstract_Particles_Positions

    type, extends(Abstract_Particles_Positions), public :: Concrete_Particles_Positions

    end type Concrete_Particles_Positions

    type, extends(Abstract_Particles_Positions), public :: Null_Particles_Positions
    contains
        procedure :: construct => Null_Particles_Positions_construct
        procedure :: destroy => Null_Particles_Positions_destroy
        procedure :: set => Null_Particles_Positions_set
        procedure :: get_num => Null_Particles_Positions_get_num
        procedure :: get => Null_Particles_Positions_get
        procedure :: add => Null_Particles_Positions_add
        procedure :: remove => Null_Particles_Positions_remove
    end type Null_Particles_Positions

contains

!implementation Abstract_Particles_Positions

    subroutine Abstract_Particles_Positions_construct(this, periodic_box, number)
        class(Abstract_Particles_Positions), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Particles_Number), target, intent(in) :: number

        this%periodic_box => periodic_box
        this%number => number
        if (this%number%get() == 0) then
            allocate(this%positions(num_dimensions, 1))
        else
            allocate(this%positions(num_dimensions, this%number%get()))
        end if
        this%positions = 0._DP
    end subroutine Abstract_Particles_Positions_construct

    subroutine Abstract_Particles_Positions_destroy(this)
        class(Abstract_Particles_Positions), intent(inout) :: this

        if (allocated(this%positions)) deallocate(this%positions)
        this%number => null()
        this%periodic_box => null()
    end subroutine Abstract_Particles_Positions_destroy

    subroutine Abstract_Particles_Positions_set(this, i_particle, position)
        class(Abstract_Particles_Positions), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: position(:)

        call check_3d_array("Abstract_Particles_Positions", "position", position)
        this%positions(:, i_particle) = this%periodic_box%folded(position)
    end subroutine Abstract_Particles_Positions_set

    pure function Abstract_Particles_Positions_get_num(this) result(num_positions)
        class(Abstract_Particles_Positions), intent(in) :: this
        integer :: num_positions

        num_positions = this%number%get()
    end function Abstract_Particles_Positions_get_num

    pure function Abstract_Particles_Positions_get(this, i_particle) result(position)
        class(Abstract_Particles_Positions), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: position(num_dimensions)

        position = this%positions(:, i_particle)
    end function Abstract_Particles_Positions_get

    subroutine Abstract_Particles_Positions_add(this, position)
        class(Abstract_Particles_Positions), intent(inout) :: this
        real(DP), intent(in) :: position(:)

        if (size(this%positions, 2) < this%number%get()) then
            call increase_coordinates_size(this%positions)
        end if
        call this%set(this%number%get(), position)
    end subroutine Abstract_Particles_Positions_add

    subroutine Abstract_Particles_Positions_remove(this, i_particle)
        class(Abstract_Particles_Positions), intent(inout) :: this
        integer, intent(in) :: i_particle

        call check_in_range("Abstract_Particles_Positions", this%number%get(), &
                            "i_particle", i_particle)
        if (i_particle < this%number%get()) then
            call this%set(i_particle, this%get(this%number%get()))
        end if
    end subroutine Abstract_Particles_Positions_remove

!end implementation Abstract_Particles_Positions

!implementation Null_Particles_Positions

    subroutine Null_Particles_Positions_construct(this, periodic_box, number)
        class(Null_Particles_Positions), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Particles_Number), target, intent(in) :: number
    end subroutine Null_Particles_Positions_construct

    subroutine Null_Particles_Positions_destroy(this)
        class(Null_Particles_Positions), intent(inout) :: this
    end subroutine Null_Particles_Positions_destroy

    subroutine Null_Particles_Positions_set(this, i_particle, position)
        class(Null_Particles_Positions), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: position(:)
    end subroutine Null_Particles_Positions_set

    pure function Null_Particles_Positions_get_num(this) result(num_positions)
        class(Null_Particles_Positions), intent(in) :: this
        integer :: num_positions
        num_positions = 0
    end function Null_Particles_Positions_get_num

    pure function Null_Particles_Positions_get(this, i_particle) result(position)
        class(Null_Particles_Positions), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: position(num_dimensions)
        position = 0._DP
    end function Null_Particles_Positions_get

    subroutine Null_Particles_Positions_add(this, position)
        class(Null_Particles_Positions), intent(inout) :: this
        real(DP), intent(in) :: position(:)
    end subroutine Null_Particles_Positions_add

    subroutine Null_Particles_Positions_remove(this, i_particle)
        class(Null_Particles_Positions), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_Particles_Positions_remove

!end implementation Null_Particles_Positions

end module class_particles_positions
