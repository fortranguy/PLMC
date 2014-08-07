module module_geometry

use module_types_micro, only: System_Geometry

implicit none
private
public geometry, set_geometry, print_geometry

    type(System_Geometry), protected :: geometry
    
contains

    subroutine set_geometry(geometry_input)
    
        character(len=*), intent(in) :: geometry_input
        
        geometry%bulk = .false.
        geometry%slab = .false.
        
        if (geometry_input == "bulk") geometry%bulk = .true.
        if (geometry_input == "slab") geometry%slab = .true.
    
    end subroutine set_geometry
    
    subroutine print_geometry(geometry_input)
    
        character(len=*), intent(in) :: geometry_input
        if (geometry_input == "bulk") write(*, *) "Geometry: Bulk"
        if (geometry_input == "slab") write(*, *) "Geometry: Slab"
        
    end subroutine print_geometry

end module module_geometry
