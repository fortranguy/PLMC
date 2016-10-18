module classes_component_coordinates_reader

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use classes_num_particles, only: Abstract_Num_Particles
use classes_component_coordinates, only: Abstract_Component_Coordinates

implicit none

private

    type, abstract, public :: Abstract_Component_Coordinates_Reader
        class(Abstract_Num_Particles), pointer :: num_particles => null()
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_read), deferred :: read
    end type Abstract_Component_Coordinates_Reader

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Component_Coordinates_Reader
            class(Abstract_Component_Coordinates_Reader), intent(inout) :: this
        end subroutine Abstract_destroy

        !> @todo integrity check?
        subroutine Abstract_read(this, coordinates_unit, num_particles)
        import :: Abstract_Component_Coordinates_Reader
            class(Abstract_Component_Coordinates_Reader), intent(in) :: this
            integer, intent(in) :: coordinates_unit
            integer, intent(in) :: num_particles
        end subroutine Abstract_read

    end interface

    type, extends(Abstract_Component_Coordinates_Reader), public :: &
        Concrete_Component_Coordinates_Reader
    private
        class(Abstract_Component_Coordinates), pointer :: positions => null()
        class(Abstract_Component_Coordinates), pointer :: orientations => null()
    contains
        procedure :: construct => Coordinates_construct
        procedure :: destroy => Coordinates_destroy
        procedure :: read => Coordinates_read
    end type Concrete_Component_Coordinates_Reader

    type, extends(Abstract_Component_Coordinates_Reader), public :: &
        Concrete_Component_Positions_Reader
    private
        class(Abstract_Component_Coordinates), pointer :: positions => null()
    contains
        procedure :: construct => Positions_construct
        procedure :: destroy => Positions_destroy
        procedure :: read => Positions_read
    end type Concrete_Component_Positions_Reader

    type, extends(Abstract_Component_Coordinates_Reader), public :: &
        Concrete_Component_Orientations_Reader
    private
        class(Abstract_Component_Coordinates), pointer :: orientations => null()
    contains
        procedure :: construct => Orientations_construct
        procedure :: destroy => Orientations_destroy
        procedure :: read => Orientations_read
    end type Concrete_Component_Orientations_Reader

    type, extends(Abstract_Component_Coordinates_Reader), public :: &
        Null_Component_Coordinates_Reader
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: read => Null_read
    end type Null_Component_Coordinates_Reader

    type, public :: Component_Coordinates_Reader_wrapper
        class(Abstract_Component_Coordinates_Reader), allocatable :: reader
    end type Component_Coordinates_Reader_wrapper

contains

!implementation Concrete_Component_Coordinates_Reader

    subroutine Coordinates_construct(this, num_particles, positions, orientations)
        class(Concrete_Component_Coordinates_Reader), intent(out) :: this
        class(Abstract_Num_Particles), target, intent(in) :: num_particles
        class(Abstract_Component_Coordinates), target, intent(in) :: positions, orientations

        this%num_particles => num_particles
        this%positions => positions
        this%orientations => orientations
    end subroutine Coordinates_construct

    subroutine Coordinates_destroy(this)
        class(Concrete_Component_Coordinates_Reader), intent(inout) :: this

        this%orientations => null()
        this%positions => null()
        this%num_particles => null()
    end subroutine Coordinates_destroy

    subroutine Coordinates_read(this, coordinates_unit, num_particles)
        class(Concrete_Component_Coordinates_Reader), intent(in) :: this
        integer, intent(in) :: coordinates_unit
        integer, intent(in) :: num_particles

        real(DP), dimension(:, :), allocatable :: positions, orientations
        integer :: i_component, i_particle

        call this%num_particles%set(num_particles)
        allocate(positions(num_dimensions, this%num_particles%get()))
        allocate(orientations(num_dimensions, this%num_particles%get()))
        do i_particle = 1, size(positions, 2)
            read(coordinates_unit, *) i_component, positions(:, i_particle), &
                orientations(:, i_particle)
        end do
        call this%positions%set_all(positions)
        call this%orientations%set_all(orientations)
    end subroutine Coordinates_read

!end implementation Concrete_Component_Coordinates_Reader

!implementation Concrete_Component_Positions_Reader

    subroutine Positions_construct(this, num_particles, positions)
        class(Concrete_Component_Positions_Reader), intent(out) :: this
        class(Abstract_Num_Particles), target, intent(in) :: num_particles
        class(Abstract_Component_Coordinates), target, intent(in) :: positions

        this%num_particles => num_particles
        this%positions => positions
    end subroutine Positions_construct

    subroutine Positions_destroy(this)
        class(Concrete_Component_Positions_Reader), intent(inout) :: this

        this%positions => null()
        this%num_particles => null()
    end subroutine Positions_destroy

    subroutine Positions_read(this, coordinates_unit, num_particles)
        class(Concrete_Component_Positions_Reader), intent(in) :: this
        integer, intent(in) :: coordinates_unit
        integer, intent(in) :: num_particles

        real(DP), dimension(:, :), allocatable :: positions
        integer :: i_component, i_particle

        call this%num_particles%set(num_particles)
        allocate(positions(num_dimensions, this%num_particles%get()))
        do i_particle = 1, size(positions, 2)
            read(coordinates_unit, *) i_component, positions(:, i_particle)
        end do
        call this%positions%set_all(positions)
    end subroutine Positions_read

!end implementation Concrete_Component_Positions_Reader

!implementation Concrete_Component_Orientations_Reader

    subroutine Orientations_construct(this, num_particles, orientations)
        class(Concrete_Component_Orientations_Reader), intent(out) :: this
        class(Abstract_Num_Particles), target, intent(in) :: num_particles
        class(Abstract_Component_Coordinates), target, intent(in) :: orientations

        this%num_particles => num_particles
        this%orientations => orientations
    end subroutine Orientations_construct

    subroutine Orientations_destroy(this)
        class(Concrete_Component_Orientations_Reader), intent(inout) :: this

        this%orientations => null()
        this%num_particles => null()
    end subroutine Orientations_destroy

    subroutine Orientations_read(this, coordinates_unit, num_particles)
        class(Concrete_Component_Orientations_Reader), intent(in) :: this
        integer, intent(in) :: coordinates_unit
        integer, intent(in) :: num_particles

        real(DP), dimension(:, :), allocatable :: orientations
        integer :: i_component, i_particle

        call this%num_particles%set(num_particles)
        allocate(orientations(num_dimensions, this%num_particles%get()))
        do i_particle = 1, size(orientations, 2)
            read(coordinates_unit, *) i_component, orientations(:, i_particle)
        end do
        call this%orientations%set_all(orientations)
    end subroutine Orientations_read

!end implementation Concrete_Component_Orientations_Reader

!implementation Null_Component_Coordinates_Reader

    subroutine Null_construct(this)
        class(Null_Component_Coordinates_Reader), intent(out) :: this
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Component_Coordinates_Reader), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_read(this, coordinates_unit, num_particles)
        class(Null_Component_Coordinates_Reader), intent(in) :: this
        integer, intent(in) :: coordinates_unit
        integer, intent(in) :: num_particles
    end subroutine Null_read

!end implementation Null_Component_Coordinates_Reader

end module classes_component_coordinates_reader
