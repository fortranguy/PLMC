module procedures_dipolar_spheres

use, intrinsic :: iso_fortran_env, only: output_unit
use data_geometry, only: num_dimensions
use class_dipolar_spheres, only: Abstract_Dipolar_Spheres

implicit none

private

public print_dipolar_spheres

contains

    subroutine print_dipolar_spheres(spheres)
        class(Abstract_Dipolar_Spheres), intent(in) :: spheres
        
        write(output_unit, *) "name = ", spheres%get_name()
        write(output_unit, *) "diameter = ", spheres%get_diameter()
        write(output_unit, *) "num = ", spheres%get_num()
        write(output_unit, *) "moment_norm = ", spheres%get_moment_norm()
        write(output_unit, *)
    end subroutine print_dipolar_spheres

end module procedures_dipolar_spheres

program test_dipolar_spheres

use class_dipolar_spheres, only: Abstract_Dipolar_Spheres, Apolar_Spheres, Dipolar_Spheres
use procedures_dipolar_spheres, only: print_dipolar_spheres
use module_data, only: test_file_exists
use json_module, only: json_file, json_initialize

implicit none

    class(Abstract_Dipolar_Spheres), allocatable :: spheres1, spheres2
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename
    
    call json_initialize()
    
    data_filename = "particles.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Apolar_Spheres :: spheres1)
    call spheres1%construct(input_data, "Spheres 1")
    allocate(Dipolar_Spheres :: spheres2)
    call spheres2%construct(input_data, "Spheres 2")
    
    call print_dipolar_spheres(spheres1)
    call print_dipolar_spheres(spheres2)
    
    call spheres2%destroy()
    deallocate(spheres2)
    call spheres1%destroy()
    deallocate(spheres1)
    
    call input_data%destroy()

end program test_dipolar_spheres
