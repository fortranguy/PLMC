module procedures_periodic_box

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_geometry, only: num_dimensions
use class_periodic_box, only: Abstract_Periodic_Box

implicit none

private

public print_periodic_box

contains

    subroutine print_periodic_box(periodic_box)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        real(DP) :: position_1(num_dimensions), position_2(num_dimensions)

        write(output_unit, *) "size =", periodic_box%get_real_size()
        position_2 = 0.9_DP * periodic_box%get_real_size()
        write(output_unit, *)
        
        position_1 = 0._DP
        position_2 = 0.9_DP * periodic_box%get_real_size()
        
        write(output_unit, *) "vector =", periodic_box%vector(position_1, position_2)
        write(output_unit, *) "distance =",  periodic_box%distance(position_1, position_2)
        write(output_unit, *)
    end subroutine print_periodic_box    

end module procedures_periodic_box

program test_periodic_box

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_geometry, only: num_dimensions
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use procedures_periodic_box, only: print_periodic_box
use module_data, only: test_file_exists, test_data_found
use json_module, only: json_file, json_initialize

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    real(DP), allocatable :: periodic_box_size(:)
    
    call json_initialize()
     
    data_filename = "periodic_box.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(XYZ_Periodic_Box :: periodic_box)
    write(output_unit, *) "XYZ"
    data_field = "Periodic Box.real size"
    call input_data%get(data_field, periodic_box_size, found)
    call test_data_found(data_field, found)
    call periodic_box%set_real_size(periodic_box_size)
    call print_periodic_box(periodic_box)
    deallocate(periodic_box)
    
    allocate(XY_Periodic_Box :: periodic_box)
    write(output_unit, *) "XY"
    data_field = "Periodic Box.real size"
    call input_data%get(data_field, periodic_box_size, found)
    call test_data_found(data_field, found)
    call periodic_box%set_real_size(periodic_box_size)
    call print_periodic_box(periodic_box)
    deallocate(periodic_box)
    
    deallocate(data_field)
    deallocate(data_filename)    
    call input_data%destroy()

end program test_periodic_box
