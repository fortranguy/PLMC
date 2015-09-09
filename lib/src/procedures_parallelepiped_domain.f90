module procedures_parallelepiped_domain

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions

implicit none

private
public point_is_inside

contains

    pure function point_is_inside(box_origin, box_size, point)
        logical :: point_is_inside
        real(DP), intent(in) :: box_origin(:), box_size(:), point(:)
        
        real(DP), dimension(num_dimensions) :: box_corner, point_from_corner
        
        box_corner = box_origin - box_size/2._DP
        point_from_corner = point - box_corner
        point_is_inside = all(0._DP <= point_from_corner .and. point_from_corner <= box_size)
    end function point_is_inside

end module procedures_parallelepiped_domain
