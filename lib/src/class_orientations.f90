module class_orientations

use, intrinsic :: iso_fortran_env, only: DP => REAL64, error_unit
use data_geometry, only: num_dimensions
use data_precisions, only: real_zero
use procedures_errors, only: warning_continue
use procedures_checks, only: check_in_range, check_3d_array, check_norm
use class_particles_number, only: Abstract_Particles_Number
use procedures_coordinates, only: increase_coordinates_size

implicit none

private

    type, abstract, public :: Abstract_Orientations
    contains
        procedure(Abstract_Orientations_construct), deferred :: construct
        procedure(Abstract_Orientations_destroy), deferred :: destroy
        procedure(Abstract_Orientations_set), deferred :: set
        procedure(Abstract_Orientations_get), deferred :: get
        procedure(Abstract_Orientations_add), deferred :: add
        procedure(Abstract_Orientations_remove), deferred :: remove
    end type Abstract_Orientations

    abstract interface

        subroutine Abstract_Orientations_construct(this, particles_number)
        import :: Abstract_Particles_Number, Abstract_Orientations
            class(Abstract_Orientations), intent(out) :: this
            class(Abstract_Particles_Number), target, intent(in) :: particles_number
        end subroutine Abstract_Orientations_construct

        subroutine Abstract_Orientations_destroy(this)
        import :: Abstract_Orientations
            class(Abstract_Orientations), intent(inout) :: this
        end subroutine Abstract_Orientations_destroy

        subroutine Abstract_Orientations_set(this, i_particle, orientation)
        import :: DP, Abstract_Orientations
            class(Abstract_Orientations), intent(inout) :: this
            integer, intent(in) :: i_particle
            real(DP), intent(in) :: orientation(:)

        end subroutine Abstract_Orientations_set

        pure function Abstract_Orientations_get(this, i_particle) result(orientation)
        import :: DP, num_dimensions, Abstract_Orientations
            class(Abstract_Orientations), intent(in) :: this
            integer, intent(in) :: i_particle
            real(DP) :: orientation(num_dimensions)

        end function Abstract_Orientations_get

        subroutine Abstract_Orientations_add(this, orientation)
        import :: DP, num_dimensions, Abstract_Orientations
            class(Abstract_Orientations), intent(inout) :: this
            real(DP), intent(in) :: orientation(:)

        end subroutine Abstract_Orientations_add

        subroutine Abstract_Orientations_remove(this, i_particle)
        import :: Abstract_Orientations
            class(Abstract_Orientations), intent(inout) :: this
            integer, intent(in) :: i_particle

        end subroutine Abstract_Orientations_remove

    end interface

    type, extends(Abstract_Orientations), public :: Null_Orientations
    contains
        procedure :: construct => Null_Orientations_construct
        procedure :: destroy => Null_Orientations_destroy

        procedure :: set => Null_Orientations_set
        procedure :: get => Null_Orientations_get
        procedure :: add => Null_Orientations_add
        procedure :: remove => Null_Orientations_remove
    end type

    type, extends(Abstract_Orientations), public :: Concrete_Orientations
    private
        real(DP), allocatable :: orientations(:, :)
        class(Abstract_Particles_Number), pointer :: particles_number
    contains
        procedure :: construct => Concrete_Orientations_construct
        procedure :: destroy => Concrete_Orientations_destroy

        procedure :: set => Concrete_Orientations_set
        procedure :: get => Concrete_Orientations_get
        procedure :: add => Concrete_Orientations_add
        procedure :: remove => Concrete_Orientations_remove
    end type Concrete_Orientations

contains

!implementation Null_Orientations

    subroutine Null_Orientations_construct(this, particles_number)
        class(Null_Orientations), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_number

    end subroutine Null_Orientations_construct

    subroutine Null_Orientations_destroy(this)
        class(Null_Orientations), intent(inout) :: this

    end subroutine Null_Orientations_destroy

    subroutine Null_Orientations_set(this, i_particle, orientation)
        class(Null_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: orientation(:)

    end subroutine Null_Orientations_set

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

!implementation Concrete_Orientations

    subroutine Concrete_Orientations_construct(this, particles_number)
        class(Concrete_Orientations), intent(out) :: this
        class(Abstract_Particles_Number), target, intent(in) :: particles_number

        integer :: i_particle

        this%particles_number => particles_number
        if (this%particles_number%get() == 0) then
            allocate(this%orientations(num_dimensions, 1))
        else
            allocate(this%orientations(num_dimensions, this%particles_number%get()))
        end if

        do i_particle = 1, this%particles_number%get()
            this%orientations(:, i_particle) = [1.0_DP, 0._DP, 0._DP]
        end do
    end subroutine Concrete_Orientations_construct

    subroutine Concrete_Orientations_destroy(this)
        class(Concrete_Orientations), intent(inout) :: this

        if (allocated(this%orientations)) deallocate(this%orientations)
        this%particles_number => null()
    end subroutine Concrete_Orientations_destroy

    subroutine Concrete_Orientations_set(this, i_particle, orientation)
        class(Concrete_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: orientation(:)

        call check_3d_array("Concrete_Orientations", "orientation", orientation)
        call check_norm("Concrete_Orientations", "orientation", orientation)
        this%orientations(:, i_particle) = orientation / norm2(orientation)
    end subroutine Concrete_Orientations_set

    pure function Concrete_Orientations_get(this, i_particle) result(orientation)
        class(Concrete_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: orientation(num_dimensions)

        orientation = this%orientations(:, i_particle)
    end function Concrete_Orientations_get

    subroutine Concrete_Orientations_add(this, orientation)
        class(Concrete_Orientations), intent(inout) :: this
        real(DP), intent(in) :: orientation(:)

        if (size(this%orientations, 2) < this%particles_number%get()) then
            call increase_coordinates_size(this%orientations)
        end if
        call this%set(this%particles_number%get(), orientation)
    end subroutine Concrete_Orientations_add

    subroutine Concrete_Orientations_remove(this, i_particle)
        class(Concrete_Orientations), intent(inout) :: this
        integer, intent(in) :: i_particle

        call check_in_range("Concrete_Orientations", this%particles_number%get(), &
                            "i_particle", i_particle)
        if (i_particle < this%particles_number%get()) then
            call this%set(i_particle, this%get(this%particles_number%get()))
        end if
    end subroutine Concrete_Orientations_remove

!end implementation Concrete_Orientations

end module class_orientations
