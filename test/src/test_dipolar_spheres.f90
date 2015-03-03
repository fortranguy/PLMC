module procedures_dipolar_spheres

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_geometry, only: num_dimensions
use class_dipolar_spheres, only: Abstract_Dipolar_Spheres

implicit none

private

public print_dipolar_spheres, set_coordinates, print_coordinates, exchange_spheres

contains

    subroutine print_dipolar_spheres(spheres)
        class(Abstract_Dipolar_Spheres), intent(in) :: spheres
        
        write(output_unit, *) "name = ", spheres%get_name()
        write(output_unit, *) "diameter = ", spheres%get_diameter()
        write(output_unit, *) "num_spheres = ", spheres%get_num_spheres()
        write(output_unit, *) "moment_norm = ", spheres%get_moment_norm()
        write(output_unit, *)
    end subroutine print_dipolar_spheres
    
    subroutine set_coordinates(spheres)
        class(Abstract_Dipolar_Spheres), intent(inout) :: spheres
        
        integer :: i_sphere
        real(DP) :: position(num_dimensions), moment(num_dimensions)
        
        moment = [1._DP, 0._DP, 0._DP] * spheres%get_moment_norm()
        do i_sphere = 1, spheres%get_num_spheres()
            position = real(i_sphere, DP)
            call spheres%set_position(i_sphere, position)
            call spheres%set_moment(i_sphere, moment)
        end do
    end subroutine set_coordinates
    
    subroutine print_coordinates(spheres)
        class(Abstract_Dipolar_Spheres), intent(in) :: spheres
        
        integer :: i_sphere
        
        do i_sphere = 1, spheres%get_num_spheres()
            write(output_unit, *) "sphere = ", i_sphere
            write(output_unit, *) "position = ", spheres%get_position(i_sphere)
            write(output_unit, *) "moment = ", spheres%get_moment(i_sphere)
        end do
    end subroutine print_coordinates
    
    subroutine exchange_spheres(spheres)
        class(Abstract_Dipolar_Spheres), intent(inout) :: spheres
        
        real(DP) :: position(num_dimensions), moment(num_dimensions)
        
        write(output_unit, *) "Add"
        write(output_unit, *) "num_spheres = ", spheres%get_num_spheres()
        position = 1.5_DP
        moment = 0.5773502691896258 * spheres%get_moment_norm()
        call spheres%add(position, moment)
        write(output_unit, *) "num_spheres = ", spheres%get_num_spheres()
        call print_coordinates(spheres)
        write(output_unit, *)
        
        write(output_unit, *) "Remove"
        write(output_unit, *) "num_spheres = ", spheres%get_num_spheres()
        call spheres%remove(1)
        write(output_unit, *) "num_spheres = ", spheres%get_num_spheres()
        call print_coordinates(spheres)
        write(output_unit, *)
        
    end subroutine exchange_spheres

end module procedures_dipolar_spheres

program test_dipolar_spheres

use, intrinsic :: iso_fortran_env, only: output_unit
use class_dipolar_spheres, only: Abstract_Dipolar_Spheres, Apolar_Spheres, Dipolar_Spheres
use procedures_dipolar_spheres, only: print_dipolar_spheres, set_coordinates, print_coordinates, &
                                      exchange_spheres
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
    
    write(output_unit, *) spheres1%get_name()
    call set_coordinates(spheres1)
    call print_coordinates(spheres1)
    call exchange_spheres(spheres1)
    
    write(output_unit, *) spheres2%get_name()
    call set_coordinates(spheres2)
    call print_coordinates(spheres2)
    call exchange_spheres(spheres2)
    
    call spheres2%destroy()
    deallocate(spheres2)
    call spheres1%destroy()
    deallocate(spheres1)
    
    call input_data%destroy()

end program test_dipolar_spheres
