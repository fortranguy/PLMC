program test_random_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use class_box_geometry, only: Abstract_Box_Geometry, Bulk_Geometry
use class_dipolar_spheres, only: Abstract_Dipolar_Spheres, Apolar_Spheres
use class_random_coordinates, only: Abstract_Random_Coordinates, Random_Coordinates
use class_random_positions, only: Abstract_Random_Positions, Bulk_Random_Positions
use class_random_moments, only: Abstract_Random_Moments, Null_Random_Moments
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
     
    data_filename = "random_coordinates.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Bulk_Geometry :: box_geometry)
    call box_geometry%set(input_data)
    allocate(Apolar_Spheres :: dipolar_spheres)
    call dipolar_spheres%construct(input_data, "Spheres 1")
    
    allocate(Bulk_Random_Positions :: rand_positions)
    call rand_positions%construct(box_geometry, dipolar_spheres, 0.5_DP)
    
    allocate(Null_Random_Moments :: rand_moments)
    call rand_moments%construct(dipolar_spheres, 1._DP)
    
    allocate(Random_Coordinates :: rand_coordinates)
    call rand_coordinates%construct(rand_positions, rand_moments)
    
    call rand_coordinates%destroy()
    deallocate(rand_coordinates)
    
    call rand_positions%destroy()
    deallocate(rand_positions)
    
    call dipolar_spheres%destroy()
    deallocate(dipolar_spheres)
    deallocate(box_geometry)
    
    call input_data%destroy()

end program test_random_coordinates
