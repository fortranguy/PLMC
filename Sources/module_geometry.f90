module module_geometry

use module_types_micro, only: System_Geometry

implicit none
private
public set_geometry

    type(System_Geometry), protected :: geometry
    
contains

    subroutine set_geometry()
    
        geometry.bulk = .true.
        geometry.slab = .false.
    
    end subroutine set_geometry

end module module_geometry
