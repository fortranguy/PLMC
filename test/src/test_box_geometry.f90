program test_box_geometry

use iso_fortran_env, only: output_unit
use module_data, only: test_file_exists
use json_module, only: json_file, json_initialize
use class_box_geometry, only: Abstract_Box_Geometry, Bulk_Geometry, Slab_Geometry

implicit none

    class(Abstract_Box_Geometry), allocatable :: box_geometry
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename
    
    call json_initialize()
    
    data_filename = "box_geometry_bulk.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Bulk_Geometry :: box_geometry)
    write(output_unit, *) "Bulk"
    call box_geometry%set(input_data)
    write(output_unit, *) "size =", box_geometry%get_size()
    write(output_unit, *) "height =", box_geometry%get_height()
    write(output_unit, *) "wave =", box_geometry%get_wave()
    write(output_unit, *)
    deallocate(box_geometry)
    
    call input_data%destroy()
    
    data_filename = "box_geometry_slab.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    
    allocate(Slab_Geometry :: box_geometry)
    write(output_unit, *) "Slab"
    call box_geometry%set(input_data)
    write(output_unit, *) "size =", box_geometry%get_size()
    write(output_unit, *) "height =", box_geometry%get_height()
    write(output_unit, *) "wave =", box_geometry%get_wave()
    write(output_unit, *)
    deallocate(box_geometry)
    
    call input_data%destroy()   

end program test_box_geometry
