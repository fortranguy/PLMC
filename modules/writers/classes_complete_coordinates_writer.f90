module classes_complete_coordinates_writer

use data_strings, only: max_line_length
use procedures_checks, only: check_positive, check_string_not_empty
use classes_number_to_string, only: Concrete_Number_to_String
use classes_box_size_writer, only: Abstract_Box_Size_Writer
use procedures_box_size_writer_factory, only: box_size_writer_destroy => destroy
use classes_component_coordinates_writer, only: Component_Coordinates_Writer_Selector
use types_component_coordinates_writer_wrapper, only: Component_Coordinates_Writer_Wrapper
use procedures_component_coordinates_writer_factory, only: component_coordinates_writer_destroy => &
    destroy

implicit none

private

    type, abstract, public :: Abstract_Complete_Coordinates_Writer
    private
        class(Abstract_Box_Size_Writer), allocatable :: box_size
        type(Component_Coordinates_Writer_Wrapper), allocatable :: components_coordinates(:)
        character(len=:), allocatable :: components_coordinates_legend
        character(len=:), allocatable :: basename
        integer :: period = 0
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: write => Abstract_write
    end type Abstract_Complete_Coordinates_Writer

contains

!implementation Abstract_Complete_Coordinates_Writer

    subroutine Abstract_construct(this, box_size, components_coordinates, &
        components_coordinates_selectors, basename, period)
        class(Abstract_Complete_Coordinates_Writer), intent(out) :: this
        class(Abstract_Box_Size_Writer), intent(in) :: box_size
        type(Component_Coordinates_Writer_Wrapper), intent(in) :: components_coordinates(:)
        type(Component_Coordinates_Writer_Selector), intent(in) :: &
            components_coordinates_selectors(:)
        character(len=*), intent(in) :: basename
        integer, intent(in) :: period

        allocate(this%box_size, source=box_size)
        allocate(this%components_coordinates, source=components_coordinates)
        this%components_coordinates_legend = "# i_component"
        if (any(components_coordinates_selectors%write_positions)) then
            this%components_coordinates_legend = this%components_coordinates_legend//&
                "    position_x    position_y    position_y"
        end if
        if (any(components_coordinates_selectors%write_orientations)) then
            this%components_coordinates_legend = this%components_coordinates_legend//&
                "    orientation_x    orientation_y    orientation_z"
        end if
        call check_string_not_empty("Abstract_Complete_Coordinates_Writer: basename", basename)
        this%basename = basename
        call check_positive("Abstract_Complete_Coordinates_Writer: construct", "period", period)
        this%period = period
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Complete_Coordinates_Writer), intent(inout) :: this

        integer :: i_component

        if (allocated(this%basename)) deallocate(this%basename)
        if (allocated(this%components_coordinates_legend)) &
            deallocate(this%components_coordinates_legend)
        if (allocated(this%components_coordinates)) then
            do i_component = size(this%components_coordinates), 1, -1
                call component_coordinates_writer_destroy(this%components_coordinates(i_component)%&
                    writer)
            end do
        end if
        call box_size_writer_destroy(this%box_size)
    end subroutine Abstract_destroy

    subroutine Abstract_write(this, i_step)
        class(Abstract_Complete_Coordinates_Writer), intent(in) :: this
        integer, intent(in) :: i_step

        integer :: nums_particles(size(this%components_coordinates)), i_component
        integer :: coordinates_unit
        type(Concrete_Number_to_String) :: string

        do i_component = 1, size(this%components_coordinates)
            nums_particles(i_component) = this%components_coordinates(i_component)%writer%get_num()
        end do

        if (mod(i_step, this%period) == 0) then
            open(newunit=coordinates_unit, recl=max_line_length, file=this%basename//"_"//string%&
                get(i_step)//".xyz", action="write")
            call this%box_size%write_new(coordinates_unit)
            write(coordinates_unit, *) "# nums_particles:", nums_particles
            write(coordinates_unit, *) this%components_coordinates_legend
            do i_component = 1, size(this%components_coordinates)
                call this%components_coordinates(i_component)%writer%write_new(coordinates_unit)
            end do
            close(coordinates_unit)
        end if
    end subroutine Abstract_write

!end implementation Abstract_Complete_Coordinates_Writer

end module classes_complete_coordinates_writer
