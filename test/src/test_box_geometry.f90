module procedures_box_geometry

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_geometry, only: num_dimensions
use class_box_geometry, only: Abstract_Box_Geometry

implicit none

private

public print_box

contains

    subroutine print_box(box)
        class(Abstract_Box_Geometry), intent(in) :: box

        real(DP) :: position1(num_dimensions), position2(num_dimensions)

        write(output_unit, *) "size =", box%get_size()
        position2 = 0.9_DP * box%get_size()
        write(output_unit, *) "height =", box%get_height()
        write(output_unit, *) "wave =", box%get_wave()
        write(output_unit, *)
        
        position1 = 0._DP
        position2(1:2) = 0.9_DP * reshape(box%get_size(), [2])
        position2(3) = 0.9_DP * box%get_height()
        
        write(output_unit, *) "vector_PBC =", box%vector_PBC(position1, position2)
        write(output_unit, *) "distance_PBC =",  box%distance_PBC(position1, position2)
        write(output_unit, *)
    end subroutine print_box    

end module procedures_box_geometry

program test_box

use, intrinsic :: iso_fortran_env, only: output_unit
use class_box_geometry, only: Abstract_Box_Geometry, Bulk_Geometry, Slab_Geometry
use procedures_box_geometry, only: print_box
use module_data, only: test_file_exists
use json_module, only: json_file, json_initialize

implicit none

    class(Abstract_Box_Geometry), allocatable :: box
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename
    
    call json_initialize()
     
    data_filename = "box_bulk.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Bulk_Geometry :: box)
    write(output_unit, *) "Bulk"
    call box%set(input_data)
    call print_box(box)
    deallocate(box)
    
    call input_data%destroy()
    
    data_filename = "box_slab.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Slab_Geometry :: box)
    write(output_unit, *) "Slab"
    call box%set(input_data)
    call print_box(box)
    deallocate(box)
    
    deallocate(data_filename)    
    call input_data%destroy()

end program test_box
