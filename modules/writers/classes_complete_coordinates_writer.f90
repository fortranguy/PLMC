module classes_complete_coordinates_writer

use data_strings, only: max_line_length
use procedures_checks, only: check_positive, check_string_not_empty
use classes_number_to_string, only: Concrete_Number_to_String
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_component_coordinates_writer, only: Component_Coordinates_Writer_Selector
use types_component_coordinates_writer_wrapper, only: Component_Coordinates_Writer_Wrapper
use procedures_component_coordinates_writer_factory, only: component_coordinates_writer_destroy => &
    destroy

implicit none

private

    type, abstract, public :: Abstract_Complete_Coordinates_Writer
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        type(Component_Coordinates_Writer_Wrapper), allocatable :: components_coordinates(:)
        character(len=:), allocatable :: components_legend
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

    subroutine Abstract_construct(this, periodic_box, components_coordinates, components_selector,&
        basename, period)
        class(Abstract_Complete_Coordinates_Writer), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        type(Component_Coordinates_Writer_Wrapper), intent(in) :: components_coordinates(:)
        type(Component_Coordinates_Writer_Selector), intent(in) :: components_selector
        character(len=*), intent(in) :: basename
        integer, intent(in) :: period

        this%periodic_box => periodic_box
        allocate(this%components_coordinates, source=components_coordinates)
        this%components_legend = "# i_component"
        if (components_selector%write_positions) then
            this%components_legend = this%components_legend//&
                "    position_x    position_y    position_z"
        end if
        if (components_selector%write_orientations) then
            this%components_legend = this%components_legend//&
                "    orientation_x    orientation_y    orientation_z"
        end if
        call check_string_not_empty("Abstract_Complete_Coordinates_Writer: basename", basename)
        this%basename = basename
        call check_positive("Abstract_Complete_Coordinates_Writer: construct", "period", period)
        this%period = period
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Complete_Coordinates_Writer), intent(inout) :: this

        if (allocated(this%basename)) deallocate(this%basename)
        if (allocated(this%components_legend)) deallocate(this%components_legend)
        call component_coordinates_writer_destroy(this%components_coordinates)
        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_write(this, i_step)
        class(Abstract_Complete_Coordinates_Writer), intent(in) :: this
        integer, intent(in) :: i_step

        integer :: nums_particles(size(this%components_coordinates)), i_component
        integer :: coordinates_unit
        type(Concrete_Number_to_String) :: string

        if (mod(i_step, this%period) /= 0) return

        do i_component = 1, size(this%components_coordinates)
            nums_particles(i_component) = this%components_coordinates(i_component)%writer%get_num()
        end do

        open(newunit=coordinates_unit, recl=max_line_length, file=this%basename//"_"//string%&
            get(i_step)//".xyz", action="write")
        write(coordinates_unit, *) "# box_size:", this%periodic_box%get_size()
        write(coordinates_unit, *) "# nums_particles:", nums_particles
        write(coordinates_unit, *) this%components_legend
        do i_component = 1, size(this%components_coordinates)
            call this%components_coordinates(i_component)%writer%write(coordinates_unit)
        end do
        close(coordinates_unit)
    end subroutine Abstract_write

!end implementation Abstract_Complete_Coordinates_Writer

!implementation Null_Complete_Coordinates_Writer

    subroutine Null_construct(this, periodic_box, components_coordinates, components_selector, &
        basename, period)
        class(Null_Complete_Coordinates_Writer), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        type(Component_Coordinates_Writer_Wrapper), intent(in) :: components_coordinates(:)
        type(Component_Coordinates_Writer_Selector), intent(in) :: components_selector
        character(len=*), intent(in) :: basename
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
