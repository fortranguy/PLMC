program test_box_geometry

use class_box_geometry, only: Abstract_Box_Geometry, Bulk_Geometry, Slab_Geometry

implicit none

    class(Abstract_Box_Geometry), allocatable :: box_geometry
    
    allocate(Bulk_Geometry :: box_geometry)
    
    deallocate(box_geometry)

end program test_box_geometry
