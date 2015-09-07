module procedures_random_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_box_geometry, only: Abstract_Box_Geometry, Bulk_Geometry, Slab_Geometry
use class_dipolar_spheres, only: Abstract_Dipolar_Spheres, Apolar_Spheres, Dipolar_Spheres
use class_random_coordinates, only: Abstract_Random_Coordinates, Random_Coordinates
use class_random_positions, only: Abstract_Random_Positions, Bulk_Random_Positions, &
                                  Slab_Random_Positions
use class_random_moments, only: Abstract_Random_Moments, Null_Random_Moments, Random_Moments
use json_module, only: json_file

implicit none

private 

public construct, destroy

contains

    subroutine construct(box, dip_spheres, rand_positions, rand_moments, &
                         rand_coordinates, input_data)
        class(Abstract_Box_Geometry), allocatable, intent(out) :: box
        class(Abstract_Dipolar_Spheres), allocatable, intent(out) :: dip_spheres
        class(Abstract_Random_Coordinates), allocatable, intent(out) :: rand_coordinates
        class(Abstract_Random_Positions), allocatable, intent(out) :: rand_positions
        class(Abstract_Random_Moments), allocatable, intent(out) :: rand_moments
        type(json_file), intent(inout) :: input_data
        
        allocate(Slab_Geometry :: box)
        call box%set(input_data)
        allocate(Dipolar_Spheres :: dip_spheres)
        call dip_spheres%construct(input_data, "Spheres 1")
        
        allocate(Slab_Random_Positions :: rand_positions)
        call rand_positions%construct(box, dip_spheres, 0.5_DP)        
        allocate(Random_Moments :: rand_moments)
        call rand_moments%construct(dip_spheres, 1._DP)        
        allocate(Random_Coordinates :: rand_coordinates)
        call rand_coordinates%construct(rand_positions, rand_moments)
    
    end subroutine construct
    
    subroutine destroy(box, dip_spheres, rand_positions, rand_moments, &
                       rand_coordinates)
        class(Abstract_Box_Geometry), allocatable, intent(inout) :: box
        class(Abstract_Dipolar_Spheres), allocatable, intent(inout) :: dip_spheres
        class(Abstract_Random_Coordinates), allocatable, intent(inout) :: rand_coordinates
        class(Abstract_Random_Positions), allocatable, intent(inout) :: rand_positions
        class(Abstract_Random_Moments), allocatable, intent(inout) :: rand_moments
    
        call rand_coordinates%destroy()
        if (allocated(rand_coordinates)) deallocate(rand_coordinates)
        
        call rand_moments%destroy()
        if (allocated(rand_moments)) deallocate(rand_moments)
        
        call rand_positions%destroy()
        if (allocated(rand_positions)) deallocate(rand_positions)
        
        call dip_spheres%destroy()
        if (allocated(dip_spheres)) deallocate(dip_spheres)
        if (allocated(box)) deallocate(box)
    
    end subroutine destroy

end module procedures_random_coordinates

program test_random_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use class_box_geometry, only: Abstract_Box_Geometry
use class_dipolar_spheres, only: Abstract_Dipolar_Spheres
use class_random_coordinates, only: Abstract_Random_Coordinates
use class_random_positions, only: Abstract_Random_Positions
use class_random_moments, only: Abstract_Random_Moments
use procedures_random_coordinates, only: construct, destroy
use module_data, only: test_file_exists
use json_module, only: json_file, json_initialize

implicit none

    class(Abstract_Box_Geometry), allocatable :: box
    class(Abstract_Dipolar_Spheres), allocatable :: dip_spheres
    class(Abstract_Random_Coordinates), allocatable :: rand_coordinates
    class(Abstract_Random_Positions), allocatable :: rand_positions
    class(Abstract_Random_Moments), allocatable :: rand_moments
    
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename
    
    call json_initialize()
     
    data_filename = "random_coordinates.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    call construct(box, dip_spheres, rand_positions, rand_moments, rand_coordinates, &
                   input_data)
                   
    write(output_unit, *) "random postion =", rand_coordinates%position()
    write(output_unit, *) "initial position =", dip_spheres%get_position(1)
    write(output_unit, *) "moved position =", rand_coordinates%move(1)
    write(output_unit, *)
    
    write(output_unit, *) "random moment =", rand_coordinates%moment()
    call dip_spheres%set_moment(1, [1._DP, 0._DP, 0._DP])
    write(output_unit, *) "moment =", dip_spheres%get_moment(1)
    write(output_unit, *) "rotated_moment =", rand_coordinates%rotation(1)
    
    call destroy(box, dip_spheres, rand_positions, rand_moments, rand_coordinates)
    
    deallocate(data_filename)
    call input_data%destroy()

end program test_random_coordinates
