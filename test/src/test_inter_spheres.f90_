program test_inter_spheres

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit

use module_data, only: test_file_exists, test_data_found
use class_sphere, only: Abstract_Sphere, Concrete_Sphere
use class_inter_spheres, only: Abstract_Inter_Spheres, Concrete_Inter_Spheres
use json_module, only: json_file, json_initialize

implicit none

    class(Abstract_Sphere), allocatable :: sphere1, sphere2
    class(Abstract_Inter_Spheres), allocatable :: inter_spheres12
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    real(DP) :: diameter, non_additivity

    call json_initialize()

    data_filename = "inter_spheres.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Concrete_Sphere :: sphere1)
    allocate(Concrete_Sphere :: sphere2)
    allocate(Concrete_Inter_Spheres :: inter_spheres12)

    data_field = "Sphere 1.diameter"
    call input_data%get(data_field, diameter, found)
    call test_data_found(data_field, found)
    call sphere1%set_diameter(diameter)
    
    data_field = "Sphere 2.diameter"
    call input_data%get(data_field, diameter, found)
    call test_data_found(data_field, found)
    call sphere2%set_diameter(diameter)
    
    data_field = "Inter Spheres 12.non additivity"
    call input_data%get(data_field, non_additivity, found)
    call test_data_found(data_field, found)
    call inter_spheres12%set_diameter(sphere1, sphere2, non_additivity)
    write(output_unit, *) "diameter =", inter_spheres12%get_diameter()

    deallocate(inter_spheres12)
    deallocate(sphere1)
    deallocate(sphere2)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_inter_spheres
