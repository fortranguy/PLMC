module class_components_energes_writer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: max_line_length
use procedures_checks, only: check_string_not_empty
use class_number_to_string, only: Abstract_Number_to_String, &
    Concrete_Number_to_String, Null_Number_to_String
use types_observables_wrapper, only: Concrete_Components_Energies

implicit none

private

    type, public :: Concrete_Components_Energies_Selector
        logical, allocatable :: with_components(:)
    end type Concrete_Components_Energies_Selector

    type :: String_Wrapper
        class(Abstract_Number_to_String), allocatable :: string
    end type String_Wrapper

    type :: Strings_Wrapper
        type(String_Wrapper), allocatable :: with_component(:)
    end type Strings_Wrapper

    type, abstract, public :: Abstract_Components_Energies_Writer
    private
        integer :: file_unit
        type(Strings_Wrapper), allocatable :: strings(:)
    contains
        procedure :: construct => Abstract_Components_Energies_Writer_construct
        procedure :: destroy => Abstract_Components_Energies_Writer_destroy
        procedure :: write => Abstract_Components_Energies_Writer_write
        procedure, private :: allocate_string => &
            Abstract_Components_Energies_Writer_allocate_strings
    end type Abstract_Components_Energies_Writer

    type, extends(Abstract_Components_Energies_Writer), public :: &
        Concrete_Components_Energies_Writer

    end type Concrete_Components_Energies_Writer

    type, extends(Abstract_Components_Energies_Writer), public :: Null_Components_Energies_Writer
    contains
        procedure :: construct => Null_Components_Energies_Writer_construct
        procedure :: destroy => Null_Components_Energies_Writer_destroy
        procedure :: write => Null_Components_Energies_Writer_write
    end type Null_Components_Energies_Writer

contains

!implementation Abstract_Components_Energies_Writer

    subroutine Abstract_Components_Energies_Writer_construct(this, filename, selector)
        class(Abstract_Components_Energies_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        type(Concrete_Components_Energies_Selector), intent(in) :: selector(:)

        character(len=:), allocatable :: legend
        integer :: file_unit !strange gfortran behaviour: otherwise writes to output_unit.

        call check_string_not_empty("Abstract_Components_Energies_Writer_construct: filename", &
            filename)
        open(newunit=file_unit, recl=max_line_length, file=filename, action="write")
        legend = "# i_step"
        call this%allocate_string(legend, selector)
        this%file_unit = file_unit
        write(this%file_unit, *) legend
    end subroutine Abstract_Components_Energies_Writer_construct

    subroutine Abstract_Components_Energies_Writer_allocate_strings(this, legend, selector)
        class(Abstract_Components_Energies_Writer), intent(out) :: this
        character(len=:), allocatable, intent(inout) :: legend
        type(Concrete_Components_Energies_Selector), intent(in) :: selector(:)

        type(Concrete_Number_to_String) :: string
        integer :: j_component, i_component

        allocate(this%strings(size(selector)))
        do j_component = 1, size(this%strings)
            allocate(this%strings(j_component)%with_component(j_component))
            do i_component = 1, size(this%strings(j_component)%with_component)
                if (selector(j_component)%with_components(i_component)) then
                    allocate(Concrete_Number_to_String :: this%strings(j_component)%&
                        with_component(i_component)%string)
                    legend = legend//"    "//string%get(i_component)//"<->"//string%get(j_component)
                else
                    allocate(Null_Number_to_String :: this%strings(j_component)%&
                        with_component(i_component)%string)
                end if
            end do
        end do
    end subroutine Abstract_Components_Energies_Writer_allocate_strings

    subroutine Abstract_Components_Energies_Writer_destroy(this)
        class(Abstract_Components_Energies_Writer), intent(inout) :: this

        integer :: i_component

        if (allocated(this%strings)) then
            do i_component = size(this%strings), 1, -1
                if (allocated(this%strings(i_component)%with_component)) then
                    deallocate(this%strings(i_component)%with_component)
                end if
            end do
            deallocate(this%strings)
        end if
        close(this%file_unit)
    end subroutine Abstract_Components_Energies_Writer_destroy

    subroutine Abstract_Components_Energies_Writer_write(this, i_step, components_energies)
        class(Abstract_Components_Energies_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Concrete_Components_Energies), intent(in) :: components_energies(:)

        character(len=:), allocatable :: energies
        integer :: j_component, i_component

        energies = ""
        do j_component = 1, size(this%strings)
            do i_component = 1, size(this%strings(j_component)%with_component)
                associate(string_ij => this%strings(j_component)%with_component(i_component)%&
                    string, energy_ij => components_energies(j_component)%&
                    with_components(i_component))
                    energies = energies//string_ij%get(energy_ij)
                end associate
            end do
        end do
        write(this%file_unit, *) i_step, energies
    end subroutine Abstract_Components_Energies_Writer_write

!end implementation Abstract_Components_Energies_Writer

!implementation Null_Components_Energies_Writer

    subroutine Null_Components_Energies_Writer_construct(this, filename, selector)
        class(Null_Components_Energies_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        type(Concrete_Components_Energies_Selector), intent(in) :: selector(:)
    end subroutine Null_Components_Energies_Writer_construct

    subroutine Null_Components_Energies_Writer_destroy(this)
        class(Null_Components_Energies_Writer), intent(inout) :: this
    end subroutine Null_Components_Energies_Writer_destroy

    subroutine Null_Components_Energies_Writer_write(this, i_step, components_energies)
        class(Null_Components_Energies_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Concrete_Components_Energies), intent(in) :: components_energies(:)
    end subroutine Null_Components_Energies_Writer_write

!end implementation Null_Components_Energies_Writer

end module class_components_energes_writer
