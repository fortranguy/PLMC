module classes_component_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64, error_unit
use data_constants, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_in_range, check_array_size, check_positive, check_norm
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_component_number, only: Abstract_Component_Number
use classes_coordinates, only: Abstract_Coordinates
use procedures_coordinates_micro, only: increase_coordinates_size

implicit none

private

    type, extends(Abstract_Coordinates), abstract, public :: Abstract_Component_Coordinates
    private
        class(Abstract_Component_Number), pointer :: number => null()
        real(DP), allocatable :: coordinates(:, :)
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure :: set_all => Abstract_set_all
        !procedure :: rescale_all
        procedure(Abstract_set), deferred :: set
        procedure :: get_num => Abstract_get_num
        procedure :: get => Abstract_get
        procedure :: add => Abstract_add
        procedure :: remove => Abstract_remove
        procedure, private :: allocate_coordinates => Abstract_allocate_coordinates
    end type Abstract_Component_Coordinates

    abstract interface

        subroutine Abstract_destroy(this)
        import Abstract_Component_Coordinates
            class(Abstract_Component_Coordinates), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_set(this, i_particle, vector)
        import :: DP, Abstract_Component_Coordinates
            class(Abstract_Component_Coordinates), intent(inout) :: this
            integer, intent(in) :: i_particle
            real(DP), intent(in) :: vector(:)
        end subroutine Abstract_set

    end interface

    type, extends(Abstract_Component_Coordinates), public :: Concrete_Component_Positions
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
    contains
        procedure :: construct => Positions_construct
        procedure :: destroy => Positions_destroy
        procedure :: set => Positions_set
    end type Concrete_Component_Positions

    type, extends(Abstract_Component_Coordinates), public :: Concrete_Component_Orientations
    contains
        procedure :: construct => Orientations_construct
        procedure :: destroy => Orientations_destroy
        procedure :: set => Orientations_set
    end type Concrete_Component_Orientations

    type, extends(Abstract_Component_Coordinates), public :: Null_Component_Coordinates
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: set_all => Null_set_all
        procedure :: set => Null_set
        procedure :: get_num => Null_get_num
        procedure :: get => Null_get
        procedure :: add => Null_add
        procedure :: remove => Null_remove
    end type Null_Component_Coordinates

contains

!implementation Abstract_Component_Coordinates

    subroutine Abstract_set_all(this, coordinates)
        class(Abstract_Component_Coordinates), intent(inout) :: this
        real(DP), intent(in) :: coordinates(:, :)

        integer :: i_particle

        if (this%get_num() /= size(coordinates, 2)) then
            call error_exit("Abstract_Component_Coordinates: set_all: numbers do not match.")
        end if
        if (size(this%coordinates, 2) < this%get_num()) then
            deallocate(this%coordinates)
            allocate(this%coordinates(num_dimensions, this%get_num()))
        end if

        do i_particle = 1, this%get_num()
            call this%set(i_particle, coordinates(:, i_particle))
        end do
    end subroutine Abstract_set_all

    subroutine Abstract_allocate_coordinates(this)
        class(Abstract_Component_Coordinates), intent(inout) :: this

        if (this%number%get() == 0) then
            allocate(this%coordinates(num_dimensions, 1))
        else
            allocate(this%coordinates(num_dimensions, this%number%get()))
        end if
    end subroutine Abstract_allocate_coordinates

    pure function Abstract_get_num(this) result(num_coordinates)
        class(Abstract_Component_Coordinates), intent(in) :: this
        integer :: num_coordinates

        num_coordinates = this%number%get()
    end function Abstract_get_num

    pure function Abstract_get(this, i_particle) result(vector)
        class(Abstract_Component_Coordinates), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: vector(num_dimensions)

        vector = this%coordinates(:, i_particle)
    end function Abstract_get

    subroutine Abstract_add(this, vector)
        class(Abstract_Component_Coordinates), intent(inout) :: this
        real(DP), intent(in) :: vector(:)

        if (size(this%coordinates, 2) < this%number%get()) then
            call increase_coordinates_size(this%coordinates)
        end if
        call this%set(this%number%get(), vector)
    end subroutine Abstract_add

    subroutine Abstract_remove(this, i_particle)
        class(Abstract_Component_Coordinates), intent(inout) :: this
        integer, intent(in) :: i_particle

        call check_in_range("Abstract_remove", this%number%get(), "i_particle", i_particle)
        if (i_particle < this%number%get()) then
            call this%set(i_particle, this%get(this%number%get()))
        end if
    end subroutine Abstract_remove

!end implementation Abstract_Component_Coordinates

!implementation Concrete_Component_Positions

    subroutine Positions_construct(this, periodic_box, number)
        class(Concrete_Component_Positions), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Number), target, intent(in) :: number

        this%periodic_box => periodic_box
        this%number => number
        call this%allocate_coordinates()
    end subroutine Positions_construct

    subroutine Positions_destroy(this)
        class(Concrete_Component_Positions), intent(inout) :: this

        if (allocated(this%coordinates)) deallocate(this%coordinates)
        this%number => null()
        this%periodic_box => null()
    end subroutine Positions_destroy

    subroutine Positions_set(this, i_particle, vector)
        class(Concrete_Component_Positions), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: vector(:)

        call check_in_range("Concrete_Component_Positions: set", this%number%get(), "i_particle", &
            i_particle)
        call check_array_size("Concrete_Component_Positions: set", "vector", vector, num_dimensions)
        this%coordinates(:, i_particle) = this%periodic_box%folded(vector)
    end subroutine Positions_set

!end implementation Concrete_Component_Positions

!implementation Concrete_Component_Orientations

    subroutine Orientations_construct(this, number)
        class(Concrete_Component_Orientations), intent(out) :: this
        class(Abstract_Component_Number), target, intent(in) :: number

        integer :: i_particle

        this%number => number
        call this%allocate_coordinates()
    end subroutine Orientations_construct

    subroutine Orientations_destroy(this)
        class(Concrete_Component_Orientations), intent(inout) :: this

        if (allocated(this%coordinates)) deallocate(this%coordinates)
        this%number => null()
    end subroutine Orientations_destroy

    subroutine Orientations_set(this, i_particle, vector)
        class(Concrete_Component_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: vector(:)

        call check_in_range("Concrete_Component_Orientations: set", this%number%get(), &
            "i_particle", i_particle)
        call check_array_size("Concrete_Component_Orientations: set", "vector", vector, &
            num_dimensions)
        call check_positive("Concrete_Component_Orientations: set", "norm2(vector)", norm2(vector))
        call check_norm("Concrete_Component_Orientations: set", "vector", vector)
        this%coordinates(:, i_particle) = vector / norm2(vector)
    end subroutine Orientations_set

!end implementation Concrete_Component_Orientations

!implementation Null_Component_Coordinates

    subroutine Null_construct(this)
        class(Null_Component_Coordinates), intent(out) :: this
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Component_Coordinates), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_set_all(this, coordinates)
        class(Null_Component_Coordinates), intent(inout) :: this
        real(DP), intent(in) :: coordinates(:, :)
    end subroutine Null_set_all

    subroutine Null_set(this, i_particle, vector)
        class(Null_Component_Coordinates), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: vector(:)
    end subroutine Null_set

    pure function Null_get_num(this) result(num_coordinates)
        class(Null_Component_Coordinates), intent(in) :: this
        integer :: num_coordinates
        num_coordinates = 0
    end function Null_get_num

    pure function Null_get(this, i_particle) result(position)
        class(Null_Component_Coordinates), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: position(num_dimensions)
        position = 0._DP
    end function Null_get

    subroutine Null_add(this, vector)
        class(Null_Component_Coordinates), intent(inout) :: this
        real(DP), intent(in) :: vector(:)
    end subroutine Null_add

    subroutine Null_remove(this, i_particle)
        class(Null_Component_Coordinates), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_remove

!end implementation Null_Component_Coordinates

end module classes_component_coordinates
