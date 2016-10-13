module classes_complete_coordinates_writer

use data_strings, only: max_line_length
use procedures_checks, only: check_positive, check_string_not_empty
use classes_number_to_string, only: Concrete_Number_to_String
use types_string_wrapper, only: String_Wrapper
use procedures_string_factory, only: string_destroy => destroy
use classes_periodic_box, only: Abstract_Periodic_Box
use types_component_coordinates_writer_selector, only: Component_Coordinates_Writer_Selector
use types_component_coordinates_writer_wrapper, only: Component_Coordinates_Writer_Wrapper
use procedures_component_coordinates_writer_factory, only: component_coordinates_writer_destroy => &
    destroy

implicit none

private

    type, abstract, public :: Abstract_Complete_Coordinates_Writer
    private
        class(Abstract_Periodic_Box), pointer :: periodic_boxes(:) => null()
        type(Component_Coordinates_Writer_Wrapper), allocatable :: components_coordinates(:, :)
        character(len=:), allocatable :: components_legend
        type(String_Wrapper), allocatable :: paths(:)
        character(len=:), allocatable :: basename
        integer :: period = 0
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: write => Abstract_write
    end type Abstract_Complete_Coordinates_Writer

    type, extends(Abstract_Complete_Coordinates_Writer), public :: &
        Concrete_Complete_Coordinates_Writer

    end type Concrete_Complete_Coordinates_Writer

    type, extends(Abstract_Complete_Coordinates_Writer), public :: Null_Complete_Coordinates_Writer
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: write => Null_write
    end type Null_Complete_Coordinates_Writer

contains

!implementation Abstract_Complete_Coordinates_Writer

    subroutine Abstract_construct(this, paths, basename, periodic_boxes, components_coordinates, &
        coordinates_selector, period)
        class(Abstract_Complete_Coordinates_Writer), intent(out) :: this
        type(String_Wrapper), intent(in) :: paths(:)
        character(len=*), intent(in) :: basename
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_boxes(:)
        type(Component_Coordinates_Writer_Wrapper), intent(in) :: components_coordinates(:, :)
        type(Component_Coordinates_Writer_Selector), intent(in) :: coordinates_selector
        integer, intent(in) :: period

        allocate(this%paths, source=paths)
        call check_string_not_empty("Abstract_Complete_Coordinates_Writer: basename", basename)
        this%basename = basename
        this%periodic_boxes => periodic_boxes
        allocate(this%components_coordinates, source=components_coordinates)
        this%components_legend = "# i_component"
        if (coordinates_selector%write_positions) then
            this%components_legend = this%components_legend//&
                "    position_x    position_y    position_z"
        end if
        if (coordinates_selector%write_orientations) then
            this%components_legend = this%components_legend//&
                "    orientation_x    orientation_y    orientation_z"
        end if
        call check_positive("Abstract_Complete_Coordinates_Writer: construct", "period", period)
        this%period = period
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Complete_Coordinates_Writer), intent(inout) :: this

        if (allocated(this%components_legend)) deallocate(this%components_legend)
        call component_coordinates_writer_destroy(this%components_coordinates)
        this%periodic_boxes => null()
        if (allocated(this%basename)) deallocate(this%basename)
        call string_destroy(this%paths)
    end subroutine Abstract_destroy

    subroutine Abstract_write(this, i_step)
        class(Abstract_Complete_Coordinates_Writer), intent(in) :: this
        integer, intent(in) :: i_step

        integer :: i_box, i_component
        integer :: nums_particles(size(this%components_coordinates, 1))
        integer :: coordinates_unit
        type(Concrete_Number_to_String) :: string

        if (mod(i_step, this%period) /= 0) return

        do i_box = 1, size(this%components_coordinates, 2)
            do i_component = 1, size(this%components_coordinates, 1)
                nums_particles(i_component) = this%components_coordinates(i_component, i_box)%&
                    writer%get_num()
            end do

            open(newunit=coordinates_unit, recl=max_line_length, file=this%paths(i_box)%string//&
                this%basename//"_"//string%get(i_step)//".xyz", action="write")
            write(coordinates_unit, *) "# box_size:", this%periodic_boxes(i_box)%get_size()
            write(coordinates_unit, *) "# nums_particles:", nums_particles
            write(coordinates_unit, *) this%components_legend
            do i_component = 1, size(this%components_coordinates, 1)
                call this%components_coordinates(i_component, i_box)%writer%write(coordinates_unit)
            end do
            close(coordinates_unit)
        end do
    end subroutine Abstract_write

!end implementation Abstract_Complete_Coordinates_Writer

!implementation Null_Complete_Coordinates_Writer

    subroutine Null_construct(this, paths, basename, periodic_boxes, components_coordinates, &
        coordinates_selector, period)
        class(Null_Complete_Coordinates_Writer), intent(out) :: this
        type(String_Wrapper), intent(in) :: paths(:)
        character(len=*), intent(in) :: basename
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_boxes(:)
        type(Component_Coordinates_Writer_Wrapper), intent(in) :: components_coordinates(:, :)
        type(Component_Coordinates_Writer_Selector), intent(in) :: coordinates_selector
        integer, intent(in) :: period
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Complete_Coordinates_Writer), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_write(this, i_step)
        class(Null_Complete_Coordinates_Writer), intent(in) :: this
        integer, intent(in) :: i_step
    end subroutine Null_write

!end implementation Null_Complete_Coordinates_Writer

end module classes_complete_coordinates_writer
