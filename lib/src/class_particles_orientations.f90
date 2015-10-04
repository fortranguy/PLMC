module class_particles_orientations

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, real_zero
use procedures_errors, only: warning_continue
use procedures_checks, only: check_in_range, check_3d_array, check_positive, check_norm
use class_particles_number, only: Abstract_Particles_Number
use procedures_coordinates, only: increase_coordinates_size
use class_particles_coordinates, only: Abstract_Particles_Coordinates

implicit none

private

    type, extends(Abstract_Particles_Coordinates), abstract, public :: &
        Abstract_Particles_Orientations
    private
        real(DP), allocatable :: orientations(:, :)
        class(Abstract_Particles_Number), pointer :: number
    contains
        procedure :: construct => Abstract_Particles_Orientations_construct
        procedure :: destroy => Abstract_Particles_Orientations_destroy
        procedure :: set => Abstract_Particles_Orientations_set
        procedure :: get_num => Abstract_Particles_Orientations_get_num
        procedure :: get => Abstract_Particles_Orientations_get
        procedure :: add => Abstract_Particles_Orientations_add
        procedure :: remove => Abstract_Particles_Orientations_remove
    end type Abstract_Particles_Orientations

    type, extends(Abstract_Particles_Orientations), public :: Concrete_Particles_Orientations

    end type Concrete_Particles_Orientations

    type, extends(Abstract_Particles_Orientations), public :: Null_Particles_Orientations
    contains
        procedure :: construct => Null_Particles_Orientations_construct
        procedure :: destroy => Null_Particles_Orientations_destroy
        procedure :: set => Null_Particles_Orientations_set
        procedure :: get_num => Null_Particles_Orientations_get_num
        procedure :: get => Null_Particles_Orientations_get
        procedure :: add => Null_Particles_Orientations_add
        procedure :: remove => Null_Particles_Orientations_remove
    end type

contains

!implementation Abstract_Particles_Orientations

    subroutine Abstract_Particles_Orientations_construct(this, number)
        class(Abstract_Particles_Orientations), intent(out) :: this
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
    end subroutine Abstract_Particles_Orientations_construct

    subroutine Abstract_Particles_Orientations_destroy(this)
        class(Abstract_Particles_Orientations), intent(inout) :: this

        if (allocated(this%orientations)) deallocate(this%orientations)
        this%number => null()
    end subroutine Abstract_Particles_Orientations_destroy

    subroutine Abstract_Particles_Orientations_set(this, i_particle, vector)
        class(Abstract_Particles_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: vector(:)

        call check_3d_array("Abstract_Particles_Orientations", "vector", vector)
        call check_positive("Abstract_Particles_Orientations", "norm2(vector)", norm2(vector))
        call check_norm("Abstract_Particles_Orientations", "vector", vector)
        this%orientations(:, i_particle) = vector / norm2(vector)
    end subroutine Abstract_Particles_Orientations_set

    pure function Abstract_Particles_Orientations_get_num(this) result(num_orientations)
        class(Abstract_Particles_Orientations), intent(in) :: this
        integer :: num_orientations

        num_orientations = this%number%get()
    end function Abstract_Particles_Orientations_get_num

    pure function Abstract_Particles_Orientations_get(this, i_particle) result(orientation)
        class(Abstract_Particles_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: orientation(num_dimensions)

        orientation = this%orientations(:, i_particle)
    end function Abstract_Particles_Orientations_get

    subroutine Abstract_Particles_Orientations_add(this, vector)
        class(Abstract_Particles_Orientations), intent(inout) :: this
        real(DP), intent(in) :: vector(:)

        if (size(this%orientations, 2) < this%number%get()) then
            call increase_coordinates_size(this%orientations)
        end if
        call this%set(this%number%get(), vector)
    end subroutine Abstract_Particles_Orientations_add

    subroutine Abstract_Particles_Orientations_remove(this, i_particle)
        class(Abstract_Particles_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle

        call check_in_range("Abstract_Particles_Orientations", this%number%get(), &
                            "i_particle", i_particle)
        if (i_particle < this%number%get()) then
            call this%set(i_particle, this%get(this%number%get()))
        end if
    end subroutine Abstract_Particles_Orientations_remove

!end implementation Abstract_Particles_Orientations

!implementation Null_Particles_Orientations

    subroutine Null_Particles_Orientations_construct(this, number)
        class(Null_Particles_Orientations), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: number
    end subroutine Null_Particles_Orientations_construct

    subroutine Null_Particles_Orientations_destroy(this)
        class(Null_Particles_Orientations), intent(inout) :: this
    end subroutine Null_Particles_Orientations_destroy

    subroutine Null_Particles_Orientations_set(this, i_particle, vector)
        class(Null_Particles_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: vector(:)
    end subroutine Null_Particles_Orientations_set

    pure function Null_Particles_Orientations_get_num(this) result(num_orientations)
        class(Null_Particles_Orientations), intent(in) :: this
        integer :: num_orientations
        num_orientations = 0
    end function Null_Particles_Orientations_get_num

    pure function Null_Particles_Orientations_get(this, i_particle) result(orientation)
        class(Null_Particles_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: orientation(num_dimensions)
        orientation = 0._DP
    end function Null_Particles_Orientations_get

    subroutine Null_Particles_Orientations_add(this, vector)
        class(Null_Particles_Orientations), intent(inout) :: this
        real(DP), intent(in) :: vector(:)
    end subroutine Null_Particles_Orientations_add

    subroutine Null_Particles_Orientations_remove(this, i_particle)
        class(Null_Particles_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_Particles_Orientations_remove

!end implementation Null_Particles_Orientations

end module class_particles_orientations
