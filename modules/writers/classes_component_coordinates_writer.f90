module classes_component_coordinates_writer

use data_strings, only: max_line_length
use procedures_checks, only: check_string_not_empty, check_positive
use classes_number_to_string, only: Abstract_Number_to_String, &
    Concrete_Number_to_String, Null_Number_to_String
use classes_component_coordinates, only: Abstract_Component_Coordinates
implicit none

private

    type, public :: Component_Coordinates_Writer_Selector
        logical :: write_positions = .false.
        logical :: write_orientations = .false.
    end type Component_Coordinates_Writer_Selector

    type, abstract, public :: Abstract_Component_Coordinates_Writer
    private
        integer :: i_component = 0
        class(Abstract_Component_Coordinates), pointer :: positions => null()
        class(Abstract_Number_to_String), allocatable :: string_positions
        class(Abstract_Number_to_String), allocatable :: string_orientations
        class(Abstract_Component_Coordinates), pointer :: orientations => null()
        type(Concrete_Number_to_String) :: string_step
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: get_num => Abstract_get_num
        procedure :: write => Abstract_write
    end type Abstract_Component_Coordinates_Writer

    type, extends(Abstract_Component_Coordinates_Writer), public :: &
        Concrete_Component_Coordinates_Writer

    end type Concrete_Component_Coordinates_Writer

    type, extends(Abstract_Component_Coordinates_Writer), public :: &
        Null_Component_Coordinates_Writer
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: get_num => Null_get_num
        procedure :: write => Null_write
    end type Null_Component_Coordinates_Writer

contains

!implementation Abstract_Component_Coordinates_Writer

    subroutine Abstract_construct(this, i_component, positions, orientations, coordinates_selector)
        class(Abstract_Component_Coordinates_Writer), intent(out) :: this
        integer, intent(in) :: i_component
        class(Abstract_Component_Coordinates), target, intent(in) :: positions, orientations
        type(Component_Coordinates_Writer_Selector), intent(in) :: coordinates_selector

        call check_positive("Abstract_Component_Coordinates_Writer: construct", "i_component", &
            i_component)
        this%i_component = i_component
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
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Component_Coordinates_Writer), intent(inout) :: this

        if (allocated(this%string_orientations)) deallocate(this%string_orientations)
        if (allocated(this%string_positions)) deallocate(this%string_positions)
        this%orientations => null()
        this%positions => null()
    end subroutine Abstract_destroy

    pure integer function Abstract_get_num(this) result(num_coordinates)
        class(Abstract_Component_Coordinates_Writer), intent(in) :: this

        num_coordinates = this%positions%get_num()
    end function Abstract_get_num

    subroutine Abstract_write(this, coordinates_unit)
        class(Abstract_Component_Coordinates_Writer), intent(in) :: this
        integer, intent(in) :: coordinates_unit

        integer :: i_particle

        do i_particle = 1, this%positions%get_num()
            write(coordinates_unit, *) this%i_component, &
                this%string_positions%get(this%positions%get(i_particle)), &
                this%string_orientations%get(this%orientations%get(i_particle))
        end do
    end subroutine Abstract_write

!end implementation Abstract_Component_Coordinates_Writer

!implementation Null_Component_Coordinates_Writer

    subroutine Null_construct(this, i_component, positions, orientations, coordinates_selector)
        class(Null_Component_Coordinates_Writer), intent(out) :: this
        integer, intent(in) :: i_component
        class(Abstract_Component_Coordinates), target, intent(in) :: positions, orientations
        type(Component_Coordinates_Writer_Selector), intent(in) :: coordinates_selector
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Component_Coordinates_Writer), intent(inout) :: this
    end subroutine Null_destroy

    pure integer function Null_get_num(this) result(num_coordinates)
        class(Null_Component_Coordinates_Writer), intent(in) :: this
        num_coordinates = 0
    end function Null_get_num

    subroutine Null_write(this, coordinates_unit)
        class(Null_Component_Coordinates_Writer), intent(in) :: this
        integer, intent(in) :: coordinates_unit
    end subroutine Null_write

!end implementation Null_Component_Coordinates_Writer

end module classes_component_coordinates_writer
