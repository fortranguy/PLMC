module class_component_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64, error_unit
use data_constants, only: num_dimensions
use procedures_checks, only: check_in_range, check_3d_array, check_positive, check_norm
use class_periodic_box, only: Abstract_Periodic_Box
use class_component_number, only: Abstract_Component_Number
use class_coordinates, only: Abstract_Coordinates
use procedures_coordinates_micro, only: increase_coordinates_size

implicit none

private

    type, extends(Abstract_Coordinates), abstract, public :: Abstract_Component_Coordinates
    private
        class(Abstract_Component_Number), pointer :: number
        real(DP), allocatable :: coordinates(:, :)
    contains
        procedure, private :: allocate_coordinates => &
            Abstract_Component_Coordinates_allocate_coordinates
        procedure(Abstract_Component_Coordinates_destroy), deferred :: destroy
        procedure(Abstract_Component_Coordinates_set), deferred :: set
        procedure :: get_num => Abstract_Component_Coordinates_get_num
        procedure :: get => Abstract_Component_Coordinates_get
        procedure :: add => Abstract_Component_Coordinates_add
        procedure :: remove => Abstract_Component_Coordinates_remove
    end type Abstract_Component_Coordinates

    abstract interface

        subroutine Abstract_Component_Coordinates_destroy(this)
        import Abstract_Component_Coordinates
            class(Abstract_Component_Coordinates), intent(inout) :: this
        end subroutine Abstract_Component_Coordinates_destroy

        subroutine Abstract_Component_Coordinates_set(this, i_particle, vector)
        import :: DP, Abstract_Component_Coordinates
            class(Abstract_Component_Coordinates), intent(inout) :: this
            integer, intent(in) :: i_particle
            real(DP), intent(in) :: vector(:)
        end subroutine Abstract_Component_Coordinates_set

    end interface

    type, extends(Abstract_Component_Coordinates), public :: Concrete_Component_Positions
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box
    contains
        procedure :: construct => Concrete_Component_Positions_construct
        procedure :: destroy => Concrete_Component_Positions_destroy
        procedure :: set => Concrete_Component_Positions_set
    end type Concrete_Component_Positions

    type, extends(Abstract_Component_Coordinates), public :: Concrete_Component_Orientations
    private
    contains
        procedure :: construct => Concrete_Component_Orientations_construct
        procedure :: destroy => Concrete_Component_Orientations_destroy
        procedure :: set => Concrete_Component_Orientations_set
    end type Concrete_Component_Orientations

    type, extends(Abstract_Component_Coordinates), public :: Null_Component_Coordinates
    contains
        procedure :: construct => Null_Component_Coordinates_construct
        procedure :: destroy => Null_Component_Coordinates_destroy
        procedure :: set => Null_Component_Coordinates_set
        procedure :: get_num => Null_Component_Coordinates_get_num
        procedure :: get => Null_Component_Coordinates_get
        procedure :: add => Null_Component_Coordinates_add
        procedure :: remove => Null_Component_Coordinates_remove
    end type Null_Component_Coordinates

contains

!implementation Abstract_Component_Coordinates

    subroutine Abstract_Component_Coordinates_allocate_coordinates(this)
        class(Abstract_Component_Coordinates), intent(inout) :: this

        if (this%number%get() == 0) then
            allocate(this%coordinates(num_dimensions, 1))
        else
            allocate(this%coordinates(num_dimensions, this%number%get()))
        end if
    end subroutine Abstract_Component_Coordinates_allocate_coordinates

    pure function Abstract_Component_Coordinates_get_num(this) result(num_coordinates)
        class(Abstract_Component_Coordinates), intent(in) :: this
        integer :: num_coordinates

        num_coordinates = this%number%get()
    end function Abstract_Component_Coordinates_get_num

    pure function Abstract_Component_Coordinates_get(this, i_particle) result(coordinate)
        class(Abstract_Component_Coordinates), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: coordinate(num_dimensions)

        coordinate = this%coordinates(:, i_particle)
    end function Abstract_Component_Coordinates_get

    subroutine Abstract_Component_Coordinates_add(this, vector)
        class(Abstract_Component_Coordinates), intent(inout) :: this
        real(DP), intent(in) :: vector(:)

        if (size(this%coordinates, 2) < this%number%get()) then
            call increase_coordinates_size(this%coordinates)
        end if
        call this%set(this%number%get(), vector)
    end subroutine Abstract_Component_Coordinates_add

    subroutine Abstract_Component_Coordinates_remove(this, i_particle)
        class(Abstract_Component_Coordinates), intent(inout) :: this
        integer, intent(in) :: i_particle

        call check_in_range("Abstract_Component_Coordinates_remove", this%number%get(), &
            "i_particle", i_particle)
        if (i_particle < this%number%get()) then
            call this%set(i_particle, this%get(this%number%get()))
        end if
    end subroutine Abstract_Component_Coordinates_remove

!end implementation Abstract_Component_Coordinates

!implementation Concrete_Component_Positions

    subroutine Concrete_Component_Positions_construct(this, periodic_box, number)
        class(Concrete_Component_Positions), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Number), target, intent(in) :: number

        this%periodic_box => periodic_box
        this%number => number
        call this%allocate_coordinates()
    end subroutine Concrete_Component_Positions_construct

    subroutine Concrete_Component_Positions_destroy(this)
        class(Concrete_Component_Positions), intent(inout) :: this

        if (allocated(this%coordinates)) deallocate(this%coordinates)
        this%number => null()
        this%periodic_box => null()
    end subroutine Concrete_Component_Positions_destroy

    subroutine Concrete_Component_Positions_set(this, i_particle, vector)
        class(Concrete_Component_Positions), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: vector(:)

        call check_in_range("Concrete_Component_Positions_set", this%number%get(), &
            "i_particle", i_particle)
        call check_3d_array("Concrete_Component_Positions_set", "vector", vector)
        this%coordinates(:, i_particle) = this%periodic_box%folded(vector)
    end subroutine Concrete_Component_Positions_set

!end implementation Concrete_Component_Positions

!implementation Concrete_Component_Orientations

    subroutine Concrete_Component_Orientations_construct(this, number)
        class(Concrete_Component_Orientations), intent(out) :: this
        class(Abstract_Component_Number), target, intent(in) :: number

        integer :: i_particle

        this%number => number
        call this%allocate_coordinates()
    end subroutine Concrete_Component_Orientations_construct

    subroutine Concrete_Component_Orientations_destroy(this)
        class(Concrete_Component_Orientations), intent(inout) :: this

        if (allocated(this%coordinates)) deallocate(this%coordinates)
        this%number => null()
    end subroutine Concrete_Component_Orientations_destroy

    subroutine Concrete_Component_Orientations_set(this, i_particle, vector)
        class(Concrete_Component_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: vector(:)

        call check_in_range("Concrete_Component_Orientations_set", this%number%get(), &
            "i_particle", i_particle)
        call check_3d_array("Concrete_Component_Orientations_set", "vector", vector)
        call check_positive("Concrete_Component_Orientations_set", "norm2(vector)", norm2(vector))
        call check_norm("Concrete_Component_Orientations_set", "vector", vector)
        this%coordinates(:, i_particle) = vector / norm2(vector)
    end subroutine Concrete_Component_Orientations_set

!end implementation Concrete_Component_Orientations

!implementation Null_Component_Coordinates

    subroutine Null_Component_Coordinates_construct(this)
        class(Null_Component_Coordinates), intent(out) :: this
    end subroutine Null_Component_Coordinates_construct

    subroutine Null_Component_Coordinates_destroy(this)
        class(Null_Component_Coordinates), intent(inout) :: this
    end subroutine Null_Component_Coordinates_destroy

    subroutine Null_Component_Coordinates_set(this, i_particle, vector)
        class(Null_Component_Coordinates), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: vector(:)
    end subroutine Null_Component_Coordinates_set

    pure function Null_Component_Coordinates_get_num(this) result(num_coordinates)
        class(Null_Component_Coordinates), intent(in) :: this
        integer :: num_coordinates
        num_coordinates = 0
    end function Null_Component_Coordinates_get_num

    pure function Null_Component_Coordinates_get(this, i_particle) result(position)
        class(Null_Component_Coordinates), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: position(num_dimensions)
        position = 0._DP
    end function Null_Component_Coordinates_get

    subroutine Null_Component_Coordinates_add(this, vector)
        class(Null_Component_Coordinates), intent(inout) :: this
        real(DP), intent(in) :: vector(:)
    end subroutine Null_Component_Coordinates_add

    subroutine Null_Component_Coordinates_remove(this, i_particle)
        class(Null_Component_Coordinates), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_Component_Coordinates_remove

!end implementation Null_Component_Coordinates

end module class_component_coordinates
