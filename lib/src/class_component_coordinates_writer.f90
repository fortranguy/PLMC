module class_component_coordinates_writer

use data_constants, only: max_line_length
use procedures_checks, only: check_string_not_empty, check_positive
use class_number_to_string, only: Abstract_Number_to_String, &
    Concrete_Number_to_String, Null_Number_to_String
use class_component_coordinates, only: Abstract_Component_Coordinates
implicit none

private

    type, public :: Concrete_Coordinates_Writer_Selector
        integer :: period
        logical :: write_positions
        logical :: write_orientations
    end type Concrete_Coordinates_Writer_Selector

    type, abstract, public :: Abstract_Coordinates_Writer
    private
        character(len=:), allocatable :: basename
        character(len=:), allocatable :: legend
        integer :: period
        type(Concrete_Number_to_String) :: string_step
        class(Abstract_Component_Coordinates), pointer :: positions => null()
        type(Concrete_Number_to_String) :: string_positions
        class(Abstract_Number_to_String), allocatable :: string_orientations
        class(Abstract_Component_Coordinates), pointer :: orientations => null()
    contains
        procedure :: construct => Abstract_Coordinates_Writer_construct
        procedure :: destroy => Abstract_Coordinates_Writer_destroy
        procedure :: write => Abstract_Coordinates_Writer_write
    end type Abstract_Coordinates_Writer

    type, extends(Abstract_Coordinates_Writer), public :: Concrete_Coordinates_Writer

    end type Concrete_Coordinates_Writer

    type, extends(Abstract_Coordinates_Writer), public :: Null_Coordinates_Writer
    contains
        procedure :: construct => Null_Coordinates_Writer_construct
        procedure :: destroy => Null_Coordinates_Writer_destroy
        procedure :: write => Null_Coordinates_Writer_write
    end type Null_Coordinates_Writer

contains

!implementation Abstract_Coordinates_Writer

    subroutine Abstract_Coordinates_Writer_construct(this, basename, positions, orientations, &
        coordinates_selector)
        class(Abstract_Coordinates_Writer), intent(out) :: this
        character(len=*), intent(in) :: basename
        class(Abstract_Component_Coordinates), target, intent(in) :: positions, orientations
        type(Concrete_Coordinates_Writer_Selector), intent(in) :: coordinates_selector

        this%positions => positions
        this%orientations => orientations
        call check_string_not_empty("Abstract_Coordinates_Writer_construct: basename", basename)
        this%basename = basename
        this%legend = "# position_x    position_y    position_z"
        call check_positive("Abstract_Coordinates_Writer_construct", &
            "coordinates_selector%period", coordinates_selector%period)
        this%period = coordinates_selector%period
        if (coordinates_selector%write_orientations) then
            allocate(Concrete_Number_to_String :: this%string_orientations)
            this%legend = this%legend//"    orientation_x    orientation_z    orientation_z"
        else
            allocate(Null_Number_to_String :: this%string_orientations)
        end if
    end subroutine Abstract_Coordinates_Writer_construct

    subroutine Abstract_Coordinates_Writer_destroy(this)
        class(Abstract_Coordinates_Writer), intent(inout) :: this

        if (allocated(this%string_orientations)) deallocate(this%string_orientations)
        if (allocated(this%legend)) deallocate(this%legend)
        if (allocated(this%basename)) deallocate(this%basename)
        this%orientations => null()
        this%positions => null()
    end subroutine Abstract_Coordinates_Writer_destroy

    subroutine Abstract_Coordinates_Writer_write(this, i_step)
        class(Abstract_Coordinates_Writer), intent(in) :: this
        integer, intent(in) :: i_step

        integer :: unit_i, i_particle

        if (mod(i_step, this%period) == 0) then
            open(newunit=unit_i, recl=max_line_length, file=this%basename//"_"//this%string_step%&
                get(i_step)//".out", action="write")
            write(unit_i, *) this%legend
            do i_particle = 1, this%positions%get_num()
                write(unit_i, *) this%string_positions%get(this%positions%get(i_particle)), &
                    this%string_orientations%get(this%orientations%get(i_particle))
            end do
            close(unit_i)
        end if
    end subroutine Abstract_Coordinates_Writer_write

!end implementation Abstract_Coordinates_Writer

!implementation Null_Coordinates_Writer

    subroutine Null_Coordinates_Writer_construct(this, basename, positions, orientations, &
        coordinates_selector)
        class(Null_Coordinates_Writer), intent(out) :: this
        character(len=*), intent(in) :: basename
        class(Abstract_Component_Coordinates), target, intent(in) :: positions, orientations
        type(Concrete_Coordinates_Writer_Selector), intent(in) :: coordinates_selector
    end subroutine Null_Coordinates_Writer_construct

    subroutine Null_Coordinates_Writer_destroy(this)
        class(Null_Coordinates_Writer), intent(inout) :: this
    end subroutine Null_Coordinates_Writer_destroy

    subroutine Null_Coordinates_Writer_write(this, i_step)
        class(Null_Coordinates_Writer), intent(in) :: this
        integer, intent(in) :: i_step
    end subroutine Null_Coordinates_Writer_write

!end implementation Null_Coordinates_Writer

end module class_component_coordinates_writer
