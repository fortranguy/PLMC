module procedures_reciprocal_lattice

use, intrinsic :: iso_fortran_env, only: output_unit
use json_module, only: json_file
use module_data, only: test_data_found
use class_periodic_box, only: Abstract_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice

implicit none

private
public print_reciprocal_lattice

contains

    subroutine print_reciprocal_lattice(input_data, periodic_box, reciprocal_lattice)
    
        type(json_file), intent(inout) :: input_data
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(inout) :: reciprocal_lattice
    
        character(len=:), allocatable :: data_field
        logical :: found
        integer, allocatable :: num(:)
        
        data_field = "Reciprocal Lattice.reci num"
        call input_data%get(data_field, num, found)
        call test_data_found(data_field, found)
        call reciprocal_lattice%construct(periodic_box, num)
        deallocate(num)
        write(output_unit, *) "reci num =", reciprocal_lattice%get_num()
        call reciprocal_lattice%destroy()
        
        deallocate(data_field)
    end subroutine print_reciprocal_lattice

end module procedures_reciprocal_lattice

program test_reciprocal_lattice

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice, Null_Reciprocal_Lattice, &
                                    Concrete_Reciprocal_Lattice
use procedures_reciprocal_lattice, only: print_reciprocal_lattice

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    class(Abstract_Reciprocal_Lattice), allocatable :: reciprocal_lattice
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    real(DP), allocatable :: periodic_box_size(:)
    
    call json_initialize()
     
    data_filename = "reciprocal_lattice.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(XYZ_Periodic_Box :: periodic_box)
    data_field = "Periodic Box.real size"
    call input_data%get(data_field, periodic_box_size, found)
    call test_data_found(data_field, found)
    call periodic_box%set_size(periodic_box_size)
    deallocate(periodic_box_size)
    
    write(output_unit, *) "Null"
    allocate(Null_Reciprocal_Lattice :: reciprocal_lattice)
    call print_reciprocal_lattice(input_data, periodic_box, reciprocal_lattice)
    deallocate(reciprocal_lattice)
    
    write(output_unit, *) "Concrete"
    allocate(Concrete_Reciprocal_Lattice :: reciprocal_lattice)
    call print_reciprocal_lattice(input_data, periodic_box, reciprocal_lattice)
    deallocate(reciprocal_lattice)
    
    deallocate(periodic_box)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_reciprocal_lattice
