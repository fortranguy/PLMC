module procedures_reciprocal_lattice_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice, Concrete_Reciprocal_Lattice, &
    Null_Reciprocal_Lattice
use procedures_property_inquirers, only: use_reciprocal_lattice

implicit none

private
public :: reciprocal_lattice_create, reciprocal_lattice_destroy

contains

    subroutine reciprocal_lattice_create(reciprocal_lattice, periodic_box, input_data, prefix)
        class(Abstract_Reciprocal_Lattice), allocatable, intent(out) :: reciprocal_lattice
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        integer, allocatable :: numbers(:)

        if (use_reciprocal_lattice(input_data, prefix)) then
            data_field = prefix//"Reciprocal Lattice.numbers"
            call input_data%get(data_field, numbers, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Reciprocal_Lattice :: reciprocal_lattice)
        else
            allocate(Null_Reciprocal_Lattice :: reciprocal_lattice)
        end if
        call reciprocal_lattice%construct(periodic_box, numbers)
    end subroutine reciprocal_lattice_create

    subroutine reciprocal_lattice_destroy(reciprocal_lattice)
        class(Abstract_Reciprocal_Lattice), allocatable, intent(inout) :: reciprocal_lattice

        if (allocated(reciprocal_lattice)) then
            call reciprocal_lattice%destroy()
            deallocate(reciprocal_lattice)
        end if
    end subroutine reciprocal_lattice_destroy


end module procedures_reciprocal_lattice_factory
