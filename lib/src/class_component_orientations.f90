module class_component_orientations

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, real_zero
use procedures_errors, only: warning_continue
use procedures_checks, only: check_in_range, check_3d_array, check_positive, check_norm
use class_particles_number, only: Abstract_Particles_Number
use procedures_coordinates_micro, only: increase_coordinates_size
use class_component_coordinates, only: Abstract_Component_Coordinates

implicit none

private

    type, extends(Abstract_Component_Coordinates), abstract, public :: &
        Abstract_Component_Orientations
    private
        real(DP), allocatable :: orientations(:, :)
        class(Abstract_Particles_Number), pointer :: number
    contains
        procedure :: construct => Abstract_Component_Orientations_construct
        procedure :: destroy => Abstract_Component_Orientations_destroy
        procedure :: set => Abstract_Component_Orientations_set
        procedure :: get_num => Abstract_Component_Orientations_get_num
        procedure :: get => Abstract_Component_Orientations_get
        procedure :: add => Abstract_Component_Orientations_add
        procedure :: remove => Abstract_Component_Orientations_remove
    end type Abstract_Component_Orientations

    type, extends(Abstract_Component_Orientations), public :: Concrete_Component_Orientations

    end type Concrete_Component_Orientations

    type, extends(Abstract_Component_Orientations), public :: Null_Component_Orientations
    contains
        procedure :: construct => Null_Component_Orientations_construct
        procedure :: destroy => Null_Component_Orientations_destroy
        procedure :: set => Null_Component_Orientations_set
        procedure :: get_num => Null_Component_Orientations_get_num
        procedure :: get => Null_Component_Orientations_get
        procedure :: add => Null_Component_Orientations_add
        procedure :: remove => Null_Component_Orientations_remove
    end type

contains

!implementation Abstract_Component_Orientations

    subroutine Abstract_Component_Orientations_construct(this, number)
        class(Abstract_Component_Orientations), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: number

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
    end subroutine Abstract_Component_Orientations_construct

    subroutine Abstract_Component_Orientations_destroy(this)
        class(Abstract_Component_Orientations), intent(inout) :: this

        if (allocated(this%orientations)) deallocate(this%orientations)
        this%number => null()
    end subroutine Abstract_Component_Orientations_destroy

    subroutine Abstract_Component_Orientations_set(this, i_particle, vector)
        class(Abstract_Component_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: vector(:)

        call check_in_range("Abstract_Component_Orientations_set", this%number%get(), &
            "i_particle", i_particle)
        call check_3d_array("Abstract_Component_Orientations_set", "vector", vector)
        call check_positive("Abstract_Component_Orientations_set", "norm2(vector)", norm2(vector))
        call check_norm("Abstract_Component_Orientations_set", "vector", vector)
        this%orientations(:, i_particle) = vector / norm2(vector)
    end subroutine Abstract_Component_Orientations_set

    pure function Abstract_Component_Orientations_get_num(this) result(num_orientations)
        class(Abstract_Component_Orientations), intent(in) :: this
        integer :: num_orientations

        num_orientations = this%number%get()
    end function Abstract_Component_Orientations_get_num

    pure function Abstract_Component_Orientations_get(this, i_particle) result(orientation)
        class(Abstract_Component_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: orientation(num_dimensions)

        orientation = this%orientations(:, i_particle)
    end function Abstract_Component_Orientations_get

    subroutine Abstract_Component_Orientations_add(this, vector)
        class(Abstract_Component_Orientations), intent(inout) :: this
        real(DP), intent(in) :: vector(:)

        if (size(this%orientations, 2) < this%number%get()) then
            call increase_coordinates_size(this%orientations)
        end if
        call this%set(this%number%get(), vector)
    end subroutine Abstract_Component_Orientations_add

    subroutine Abstract_Component_Orientations_remove(this, i_particle)
        class(Abstract_Component_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle

        call check_in_range("Abstract_Component_Orientations", this%number%get(), &
            "i_particle", i_particle)
        if (i_particle < this%number%get()) then
            call this%set(i_particle, this%get(this%number%get()))
        end if
    end subroutine Abstract_Component_Orientations_remove

!end implementation Abstract_Component_Orientations

!implementation Null_Component_Orientations

    subroutine Null_Component_Orientations_construct(this, number)
        class(Null_Component_Orientations), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: number
    end subroutine Null_Component_Orientations_construct

    subroutine Null_Component_Orientations_destroy(this)
        class(Null_Component_Orientations), intent(inout) :: this
    end subroutine Null_Component_Orientations_destroy

    subroutine Null_Component_Orientations_set(this, i_particle, vector)
        class(Null_Component_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: vector(:)
    end subroutine Null_Component_Orientations_set

    pure function Null_Component_Orientations_get_num(this) result(num_orientations)
        class(Null_Component_Orientations), intent(in) :: this
        integer :: num_orientations
        num_orientations = 0
    end function Null_Component_Orientations_get_num

    pure function Null_Component_Orientations_get(this, i_particle) result(orientation)
        class(Null_Component_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: orientation(num_dimensions)
        orientation = 0._DP
    end function Null_Component_Orientations_get

    subroutine Null_Component_Orientations_add(this, vector)
        class(Null_Component_Orientations), intent(inout) :: this
        real(DP), intent(in) :: vector(:)
    end subroutine Null_Component_Orientations_add

    subroutine Null_Component_Orientations_remove(this, i_particle)
        class(Null_Component_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_Component_Orientations_remove

!end implementation Null_Component_Orientations

end module class_component_orientations
