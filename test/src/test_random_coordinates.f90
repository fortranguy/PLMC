program test_random_coordinates

use, intrinsic :: iso_fortran_env, only: output_unit
use class_box_geometry, only: Abstract_Box_Geometry
use class_dipolar_spheres, only: Abstract_Dipolar_Spheres
use class_random_coordinates, only: Abstract_Random_Coordinates
use class_random_positions, only: Abstract_Random_Positions
use class_random_moments, only: Abstract_Random_Moments
use module_data, only: test_file_exists
use json_module, only: json_file, json_initialize

implicit none

    class(Abstract_Box_Geometry), allocatable :: box_geometry
    class(Abstract_Dipolar_Spheres), allocatable :: dipolar_spheres
    class(Abstract_Random_Coordinates), allocatable :: rand_coordinates
    class(Abstract_Random_Positions), allocatable :: rand_positions
    class(Abstract_Random_Moments), allocatable :: rand_moments
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename
    
    call json_initialize()
     
    data_filename = "box_geometry_bulk.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    call input_data%destroy()

end program test_random_coordinates
