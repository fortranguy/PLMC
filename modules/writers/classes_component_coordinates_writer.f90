module classes_component_coordinates_writer

use data_strings, only: max_line_length
use procedures_checks, only: check_string_not_empty, check_positive
use classes_number_to_string, only: Abstract_Number_to_String, &
    Concrete_Number_to_String, Null_Number_to_String
use classes_component_coordinates, only: Abstract_Component_Coordinates
implicit none

private

    type, public :: Concrete_Coordinates_Writer_Selector
        integer :: period = 0
        logical :: write_positions = .false.
        logical :: write_orientations = .false.
    end type Concrete_Coordinates_Writer_Selector

    type, abstract, public :: Abstract_Coordinates_Writer
    private
        character(len=:), allocatable :: i_component
        class(Abstract_Component_Coordinates), pointer :: positions => null()
        class(Abstract_Number_to_String), allocatable :: string_positions
        class(Abstract_Number_to_String), allocatable :: string_orientations
        class(Abstract_Component_Coordinates), pointer :: orientations => null()
        character(len=:), allocatable :: basename
        character(len=:), allocatable :: legend
        integer :: period = 0
        type(Concrete_Number_to_String) :: string_step
    contains
        procedure :: construct_new => Abstract_construct_new
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: write_new => Abstract_write_new
        procedure :: write => Abstract_write
    end type Abstract_Coordinates_Writer

    type, extends(Abstract_Coordinates_Writer), public :: Concrete_Coordinates_Writer

    end type Concrete_Coordinates_Writer

    type, extends(Abstract_Coordinates_Writer), public :: Null_Coordinates_Writer
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: write => Null_write
    end type Null_Coordinates_Writer

contains

!implementation Abstract_Coordinates_Writer

    subroutine Abstract_construct_new(this, i_component, positions, orientations, &
        coordinates_selector)
        class(Abstract_Coordinates_Writer), intent(out) :: this
        integer, intent(in) :: i_component
        class(Abstract_Component_Coordinates), target, intent(in) :: positions, orientations
        type(Concrete_Coordinates_Writer_Selector), intent(in) :: coordinates_selector

        type(Concrete_Number_to_String) :: string

        call check_positive("Abstract_Coordinates_Writer: construct", "i_component", i_component)
        this%i_component = string%get(i_component)
        this%positions => positions
        this%orientations => orientations
        if (coordinates_selector%write_positions) then
            allocate(Concrete_Number_to_String :: this%string_positions)
        else
            allocate(Null_Number_to_String :: this%string_positions)
        end if
        if (coordinates_selector%write_orientations) then
            allocate(Concrete_Number_to_String :: this%string_orientations)
        else
            allocate(Null_Number_to_String :: this%string_orientations)
        end if
    end subroutine Abstract_construct_new

    subroutine Abstract_construct(this, positions, orientations, coordinates_selector, basename)
        class(Abstract_Coordinates_Writer), intent(out) :: this
        class(Abstract_Component_Coordinates), target, intent(in) :: positions, orientations
        type(Concrete_Coordinates_Writer_Selector), intent(in) :: coordinates_selector
        character(len=*), intent(in) :: basename

        call check_positive("Abstract_Coordinates_Writer: construct", &
            "coordinates_selector%period", coordinates_selector%period)
        this%period = coordinates_selector%period
        call check_string_not_empty("Abstract_Coordinates_Writer: construct: basename", basename)
        this%basename = basename
        this%legend = "#"
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Coordinates_Writer), intent(inout) :: this

        if (allocated(this%string_orientations)) deallocate(this%string_orientations)
        if (allocated(this%string_positions)) deallocate(this%string_positions)
        if (allocated(this%legend)) deallocate(this%legend)
        if (allocated(this%basename)) deallocate(this%basename)
        this%orientations => null()
        this%positions => null()
        if (allocated(this%i_component)) deallocate(this%i_component)
    end subroutine Abstract_destroy

    subroutine Abstract_write_new(this, coordinates_unit)
        class(Abstract_Coordinates_Writer), intent(in) :: this
        integer, intent(in) :: coordinates_unit

        integer :: i_particle

        do i_particle = 1, this%positions%get_num()
            write(coordinates_unit, *) this%i_component, &
                this%string_positions%get(this%positions%get(i_particle)), &
                this%string_orientations%get(this%orientations%get(i_particle))
        end do
    end subroutine Abstract_write_new

    subroutine Abstract_write(this, i_step)
        class(Abstract_Coordinates_Writer), intent(in) :: this
        integer, intent(in) :: i_step

        integer :: unit_i

        if (mod(i_step, this%period) == 0) then
            open(newunit=unit_i, recl=max_line_length, file=this%basename//"_"//this%string_step%&
                get(i_step)//".out", action="write")
            write(unit_i, *) "# number", this%positions%get_num()
            write(unit_i, *) this%legend
            close(unit_i)
        end if
    end subroutine Abstract_write

!end implementation Abstract_Coordinates_Writer

!implementation Null_Coordinates_Writer

    subroutine Null_construct_new(this, i_component, positions, orientations, coordinates_selector)
        class(Null_Coordinates_Writer), intent(out) :: this
        integer, intent(in) :: i_component
        class(Abstract_Component_Coordinates), target, intent(in) :: positions, orientations
        type(Concrete_Coordinates_Writer_Selector), intent(in) :: coordinates_selector
    end subroutine Null_construct_new

    subroutine Null_construct(this, positions, orientations, coordinates_selector, basename)
        class(Null_Coordinates_Writer), intent(out) :: this
        class(Abstract_Component_Coordinates), target, intent(in) :: positions, orientations
        type(Concrete_Coordinates_Writer_Selector), intent(in) :: coordinates_selector
        character(len=*), intent(in) :: basename
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Coordinates_Writer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_write(this, i_step)
        class(Null_Coordinates_Writer), intent(in) :: this
        integer, intent(in) :: i_step
    end subroutine Null_write

!end implementation Null_Coordinates_Writer

end module classes_component_coordinates_writer
