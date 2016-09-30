module procedures_reciprocal_lattices_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice, Concrete_Reciprocal_Lattice, &
    Null_Reciprocal_Lattice
use procedures_environment_inquirers, only: use_reciprocal_lattice

implicit none

private
public :: create, destroy

contains

    subroutine create(reciprocal_lattices, periodic_boxes, generating_data, prefix)
        class(Abstract_Reciprocal_Lattice), allocatable, intent(out) :: reciprocal_lattices(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: i_box
        character(len=:), allocatable :: data_field
        logical :: data_found
        integer, allocatable :: numbers(:)

        if (use_reciprocal_lattice(generating_data, prefix)) then
            data_field = prefix//"Reciprocal Lattice.numbers"
            call generating_data%get(data_field, numbers, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Reciprocal_Lattice :: reciprocal_lattices(size(periodic_boxes)))
        else
            allocate(Null_Reciprocal_Lattice :: reciprocal_lattices(size(periodic_boxes)))
        end if
        do i_box = 1, size(reciprocal_lattices)
            call reciprocal_lattices(i_box)%construct(periodic_boxes(i_box), numbers)
        end do
    end subroutine create

    subroutine destroy(reciprocal_lattices)
        class(Abstract_Reciprocal_Lattice), allocatable, intent(inout) :: reciprocal_lattices(:)

        integer :: i_box

        if (allocated(reciprocal_lattices)) then
            do i_box = size(reciprocal_lattices), 1, -1
                call reciprocal_lattices(i_box)%destroy()
            end do
            deallocate(reciprocal_lattices)
        end if
    end subroutine destroy


end module procedures_reciprocal_lattices_factory
