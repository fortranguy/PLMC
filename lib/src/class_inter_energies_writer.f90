module class_inter_energes_writer

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_components, max_line_length
use procedures_checks, only: check_string_not_empty
use class_number_to_string, only: Abstract_Number_to_String, &
    Concrete_Number_to_String, Null_Number_to_String
use types_observables_wrapper, only: Concrete_Inter_Energies

implicit none

private

    type, public :: Concrete_Inter_Energy_Selector
        logical, allocatable :: with_components(:)
    end type Concrete_Inter_Energy_Selector

    type :: String_Wrapper
        class(Abstract_Number_to_String), allocatable :: string
    end type String_Wrapper

    type :: Strings_Wrapper
        type(String_Wrapper), allocatable :: with_component(:)
    end type Strings_Wrapper

    type, abstract, public :: Abstract_Inter_Energy_Writer
    private
        integer :: unit
        type(Strings_Wrapper), allocatable :: strings(:)
    contains
        procedure :: construct => Abstract_Inter_Energy_Writer_construct
        procedure :: destroy => Abstract_Inter_Energy_Writer_destroy
        procedure :: write => Abstract_Inter_Energy_Writer_write
        procedure, private :: allocate_string => Abstract_Inter_Energy_Writer_allocate_strings
    end type Abstract_Inter_Energy_Writer

    type, extends(Abstract_Inter_Energy_Writer), public :: Concrete_Inter_Energy_Writer

    end type Concrete_Inter_Energy_Writer

    type, extends(Abstract_Inter_Energy_Writer), public :: Null_Inter_Energy_Writer
    contains
        procedure :: construct => Null_Inter_Energy_Writer_construct
        procedure :: destroy => Null_Inter_Energy_Writer_destroy
        procedure :: write => Null_Inter_Energy_Writer_write
    end type Null_Inter_Energy_Writer

contains

!implementation Abstract_Inter_Energy_Writer

    subroutine Abstract_Inter_Energy_Writer_construct(this, filename, selector)
        class(Abstract_Inter_Energy_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        type(Concrete_Inter_Energy_Selector), intent(in) :: selector(:)

        character(len=:), allocatable :: legend

        call check_string_not_empty("Abstract_Inter_Energy_Writer_construct: filename", &
            filename)
        open(newunit=this%unit, recl=max_line_length, file=filename, action="write")
        legend = "# i_step"
        call this%allocate_string(legend, selector)
        write(this%unit, *) legend
        deallocate(legend)
    end subroutine Abstract_Inter_Energy_Writer_construct

    subroutine Abstract_Inter_Energy_Writer_allocate_strings(this, legend, selector)
        class(Abstract_Inter_Energy_Writer), intent(out) :: this
        character(len=:), allocatable, intent(inout) :: legend
        type(Concrete_Inter_Energy_Selector), intent(in) :: selector(:)

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
    end subroutine Abstract_Inter_Energy_Writer_allocate_strings

    subroutine Abstract_Inter_Energy_Writer_destroy(this)
        class(Abstract_Inter_Energy_Writer), intent(inout) :: this

        integer :: i_component

        if (allocated(this%strings)) then
            do i_component = size(this%strings), 1, -1
                if (allocated(this%strings(i_component)%with_component)) then
                    deallocate(this%strings(i_component)%with_component)
                end if
            end do
        end if
        close(this%unit)
    end subroutine Abstract_Inter_Energy_Writer_destroy

    subroutine Abstract_Inter_Energy_Writer_write(this, i_step, inter_energies)
        class(Abstract_Inter_Energy_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Concrete_Inter_Energies), intent(in) :: inter_energies(:)

        character(len=:), allocatable :: energies
        integer :: j_component, i_component

        energies = ""
        do j_component = 1, size(this%strings)
            do i_component = 1, size(this%strings(j_component)%with_component)
                associate(string_ij => this%strings(j_component)%with_component(i_component)%&
                    string, energy_ij => inter_energies(j_component)%with_components(i_component))
                    energies = energies//string_ij%get(energy_ij)
                end associate
            end do
        end do
        write(this%unit, *) i_step, energies
        deallocate(energies)
    end subroutine Abstract_Inter_Energy_Writer_write

!end implementation Abstract_Inter_Energy_Writer

!implementation Null_Inter_Energy_Writer

    subroutine Null_Inter_Energy_Writer_construct(this, filename, selector)
        class(Null_Inter_Energy_Writer), intent(out) :: this
        character(len=*), intent(in) :: filename
        type(Concrete_Inter_Energy_Selector), intent(in) :: selector(:)
    end subroutine Null_Inter_Energy_Writer_construct

    subroutine Null_Inter_Energy_Writer_destroy(this)
        class(Null_Inter_Energy_Writer), intent(inout) :: this
    end subroutine Null_Inter_Energy_Writer_destroy

    subroutine Null_Inter_Energy_Writer_write(this, i_step, inter_energies)
        class(Null_Inter_Energy_Writer), intent(in) :: this
        integer, intent(in) :: i_step
        type(Concrete_Inter_Energies), intent(in) :: inter_energies(:)
    end subroutine Null_Inter_Energy_Writer_write

!end implementation Null_Inter_Energy_Writer

end module class_inter_energes_writer
