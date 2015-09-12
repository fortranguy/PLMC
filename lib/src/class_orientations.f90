module class_orientations

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use data_precisions, only: real_zero
use procedures_errors, only: warning_continue
use procedures_checks, only: check_in_range, check_3d_array, check_positive, check_norm
use class_number, only: Abstract_Number
use procedures_coordinates, only: increase_coordinates_size

implicit none

private

    type, abstract, public :: Abstract_Orientations
    private
        real(DP), allocatable :: orientations(:, :)
        class(Abstract_Number), pointer :: number
    contains
        procedure :: construct => Abstract_Orientations_construct
        procedure :: destroy => Abstract_Orientations_destroy
        procedure :: set => Abstract_Orientations_set
        procedure :: get_num => Abstract_Orientations_get_num
        procedure :: get => Abstract_Orientations_get
        procedure :: add => Abstract_Orientations_add
        procedure :: remove => Abstract_Orientations_remove
    end type Abstract_Orientations

    type, extends(Abstract_Orientations), public :: Concrete_Orientations

    end type Concrete_Orientations

    type, extends(Abstract_Orientations), public :: Null_Orientations
    contains
        procedure :: construct => Null_Orientations_construct
        procedure :: destroy => Null_Orientations_destroy
        procedure :: set => Null_Orientations_set
        procedure :: get_num => Null_Orientations_get_num
        procedure :: get => Null_Orientations_get
        procedure :: add => Null_Orientations_add
        procedure :: remove => Null_Orientations_remove
    end type

contains

!implementation Abstract_Orientations

    subroutine Abstract_Orientations_construct(this, number)
        class(Abstract_Orientations), intent(out) :: this
        class(Abstract_Number), target, intent(in) :: number

        integer :: i_particle

        this%number => number
        if (this%number%get() == 0) then
            allocate(this%orientations(num_dimensions, 1))
        else
            allocate(this%orientations(num_dimensions, this%number%get()))
        end if
        do i_particle = 1, this%number%get()
            this%orientations(:, i_particle) = 0._DP
        end do
    end subroutine Abstract_Orientations_construct

    subroutine Abstract_Orientations_destroy(this)
        class(Abstract_Orientations), intent(inout) :: this

        if (allocated(this%orientations)) deallocate(this%orientations)
        this%number => null()
    end subroutine Abstract_Orientations_destroy

    subroutine Abstract_Orientations_set(this, i_particle, orientation)
        class(Abstract_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: orientation(:)

        call check_3d_array("Abstract_Orientations", "orientation", orientation)
        call check_positive("Abstract_Orientations", "norm2(orientation)", norm2(orientation))
        call check_norm("Abstract_Orientations", "orientation", orientation)
        this%orientations(:, i_particle) = orientation / norm2(orientation)
    end subroutine Abstract_Orientations_set

    pure function Abstract_Orientations_get_num(this) result(num_orientations)
        class(Abstract_Orientations), intent(in) :: this
        integer :: num_orientations

        num_orientations = this%number%get()
    end function Abstract_Orientations_get_num

    pure function Abstract_Orientations_get(this, i_particle) result(orientation)
        class(Abstract_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: orientation(num_dimensions)

        orientation = this%orientations(:, i_particle)
    end function Abstract_Orientations_get

    subroutine Abstract_Orientations_add(this, orientation)
        class(Abstract_Orientations), intent(inout) :: this
        real(DP), intent(in) :: orientation(:)

        if (size(this%orientations, 2) < this%number%get()) then
            call increase_coordinates_size(this%orientations)
        end if
        call this%set(this%number%get(), orientation)
    end subroutine Abstract_Orientations_add

    subroutine Abstract_Orientations_remove(this, i_particle)
        class(Abstract_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle

        call check_in_range("Abstract_Orientations", this%number%get(), &
                            "i_particle", i_particle)
        if (i_particle < this%number%get()) then
            call this%set(i_particle, this%get(this%number%get()))
        end if
    end subroutine Abstract_Orientations_remove

!end implementation Abstract_Orientations

!implementation Null_Orientations

    subroutine Null_Orientations_construct(this, number)
        class(Null_Orientations), intent(out) :: this
        class(Abstract_Number), target, intent(in) :: number
    end subroutine Null_Orientations_construct

    subroutine Null_Orientations_destroy(this)
        class(Null_Orientations), intent(inout) :: this
    end subroutine Null_Orientations_destroy

    subroutine Null_Orientations_set(this, i_particle, orientation)
        class(Null_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: orientation(:)
    end subroutine Null_Orientations_set

    pure function Null_Orientations_get_num(this) result(num_orientations)
        class(Null_Orientations), intent(in) :: this
        integer :: num_orientations
        num_orientations = 0
    end function Null_Orientations_get_num

    pure function Null_Orientations_get(this, i_particle) result(orientation)
        class(Null_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: orientation(num_dimensions)
        orientation = 0._DP
    end function Null_Orientations_get

    subroutine Null_Orientations_add(this, orientation)
        class(Null_Orientations), intent(inout) :: this
        real(DP), intent(in) :: orientation(:)
    end subroutine Null_Orientations_add

    subroutine Null_Orientations_remove(this, i_particle)
        class(Null_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_Orientations_remove

!end implementation Null_Orientations

end module class_orientations
