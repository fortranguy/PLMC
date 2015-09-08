program test_reciprocal_lattice

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice, Null_Reciprocal_Lattice, &
                                    Concrete_Reciprocal_Lattice
use module_data, only: test_file_exists, test_data_found
use json_module, only: json_file, json_initialize

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    class(Abstract_Reciprocal_Lattice), allocatable :: reciprocal_lattice
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    real(DP), allocatable :: periodic_box_size(:)
    integer, allocatable :: reci_num(:)
    
    call json_initialize()
     
    data_filename = "reciprocal_lattice.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(XYZ_Periodic_Box :: periodic_box)
    data_field = "Periodic Box.real size"
    call input_data%get(data_field, periodic_box_size, found)
    call test_data_found(data_field, found)
    call periodic_box%set_real_size(periodic_box_size)
    deallocate(periodic_box_size)
    
    allocate(Null_Reciprocal_Lattice :: reciprocal_lattice)
    write(output_unit, *) "Null"
    data_field = "Reciprocal Lattice.reci num"
    call input_data%get(data_field, reci_num, found)
    call test_data_found(data_field, found)
    call reciprocal_lattice%construct(periodic_box, reci_num)
    deallocate(reci_num)
    write(output_unit, *) "reci num =", reciprocal_lattice%get_reci_num()
    call reciprocal_lattice%destroy()
    deallocate(reciprocal_lattice)
    
    allocate(Concrete_Reciprocal_Lattice :: reciprocal_lattice)
    write(output_unit, *) "Concrete"
    data_field = "Reciprocal Lattice.reci num"
    call input_data%get(data_field, reci_num, found)
    call test_data_found(data_field, found)
    call reciprocal_lattice%construct(periodic_box, reci_num)
    deallocate(reci_num)
    write(output_unit, *) "reci num =", reciprocal_lattice%get_reci_num()
    call reciprocal_lattice%destroy()
    deallocate(reciprocal_lattice)
    
    deallocate(periodic_box)
    deallocate(data_field)
    deallocate(data_filename)    
    call input_data%destroy()

end program test_reciprocal_lattice
