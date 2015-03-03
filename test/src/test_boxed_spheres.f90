module procedures_boxed_spheres

use class_box_geometry, only: Abstract_Box_Geometry, Bulk_Geometry, Slab_Geometry
use class_dipolar_spheres, only: Abstract_Dipolar_Spheres, Apolar_Spheres, Dipolar_Spheres
use class_boxed_spheres, only: Abstract_Boxed_Spheres, Bulk_Boxed_Spheres, Slab_Boxed_Spheres
use json_module, only: json_file

implicit none

private

public construct, destroy

contains

    subroutine construct(box, dip_spheres, boxed_spheres, input_data)
        class(Abstract_Box_Geometry), allocatable, intent(out) :: box
        class(Abstract_Dipolar_Spheres), allocatable, intent(out) :: dip_spheres
        class(Abstract_Boxed_Spheres), allocatable, intent(out) :: boxed_spheres
        type(json_file), intent(inout) :: input_data

        allocate(Bulk_Geometry :: box)
        call box%set(input_data)
        allocate(Dipolar_Spheres :: dip_spheres)
        call dip_spheres%construct(input_data, "Spheres 1")
        allocate(Bulk_Boxed_Spheres :: boxed_spheres)
        call boxed_spheres%construct(box, dip_spheres)
        
    end subroutine construct

    subroutine destroy(box, dip_spheres, boxed_spheres)
        class(Abstract_Box_Geometry), allocatable, intent(inout) :: box
        class(Abstract_Dipolar_Spheres), allocatable, intent(inout) :: dip_spheres
        class(Abstract_Boxed_Spheres), allocatable, intent(inout) :: boxed_spheres

        call boxed_spheres%destroy()

        call dip_spheres%destroy()
        if (allocated(dip_spheres)) deallocate(dip_spheres)

        if (allocated(box)) deallocate(box)
    end subroutine destroy

end module procedures_boxed_spheres

program test_boxed_spheres

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_geometry, only: num_dimensions
use class_box_geometry, only: Abstract_Box_Geometry
use class_dipolar_spheres, only: Abstract_Dipolar_Spheres
use class_boxed_spheres, only: Abstract_Boxed_Spheres
use procedures_boxed_spheres, only: construct, destroy
use module_data, only: test_file_exists
use json_module, only: json_file, json_initialize

implicit none

    class(Abstract_Box_Geometry), allocatable :: box
    class(Abstract_Dipolar_Spheres), allocatable :: dip_spheres
    class(Abstract_Boxed_Spheres), allocatable :: boxed_spheres
    real(DP) :: position(num_dimensions)
    
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()

    data_filename = "boxed_spheres.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    call construct(box, dip_spheres, boxed_spheres, input_data)

    write(output_unit, *) "accessible volume =", boxed_spheres%accessible_volume()
    position = 1.5_DP * box%get_size()
    write(output_unit, *) "position =", position ," is inside =", boxed_spheres%is_inside(position)
    call boxed_spheres%fold(position)
    write(output_unit, *) "folded position =", position

    call destroy(box, dip_spheres, boxed_spheres)

    call input_data%destroy()

end program test_boxed_spheres
