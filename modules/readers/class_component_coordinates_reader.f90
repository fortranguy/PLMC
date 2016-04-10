module class_component_coordinates_reader

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_component_number, only: Abstract_Component_Number
use class_component_coordinates, only: Abstract_Component_Coordinates
use procedures_coordinates_reader, only: create_coordinates_from_file

implicit none

private

    type, public :: Concrete_Coordinates_Reader_Selector
        logical :: read_positions
        logical :: read_orientations
    end type Concrete_Coordinates_Reader_Selector

    type, abstract, public :: Abstract_Coordinates_Reader
    private
        class(Abstract_Component_Number), pointer :: number => null()
        class(Abstract_Component_Coordinates), pointer :: positions => null()
        class(Abstract_Component_Coordinates), pointer :: orientations => null()
        logical :: read_orientations
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: read => Abstract_read
    end type Abstract_Coordinates_Reader

contains

!implementation Abstract_Coordinates_Reader

    subroutine Abstract_construct(this, number, positions, orientations, read_orientations)
        class(Abstract_Coordinates_Reader), intent(out) :: this
        class(Abstract_Component_Number), target, intent(in) :: number
        class(Abstract_Component_Coordinates), target, intent(in) :: positions
        class(Abstract_Component_Coordinates), target, intent(in) :: orientations
        logical, intent(in) :: read_orientations

        this%number => number
        this%positions => positions
        this%orientations => orientations
        this%read_orientations = read_orientations
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Coordinates_Reader), intent(inout) :: this

        this%orientations => null()
        this%positions => null()
        this%number => null()
    end subroutine Abstract_destroy

    subroutine Abstract_read(this, filename)
        class(Abstract_Coordinates_Reader), intent(in) :: this
        character(len=*), intent(in) :: filename

        real(DP), dimension(:, :), allocatable :: positions, orientations

        call create_coordinates_from_file(positions, orientations, filename, this%read_orientations)
        call this%number%set(size(positions, 2))
        call this%positions%set_all(positions)
        call this%orientations%set_all(orientations)
    end subroutine Abstract_read

!end implementation Abstract_Coordinates_Reader

end module class_component_coordinates_reader
