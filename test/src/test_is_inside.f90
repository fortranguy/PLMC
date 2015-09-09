program test_point_is_inside

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_geometry, only: num_dimensions
use procedures_parallelepiped_domain, only: point_is_inside

implicit none

    real(DP), dimension(num_dimensions):: box_origin, box_size, point
    
    box_origin = [5._DP, 5._DP, 5._DP]
    box_size = [10._DP, 10._DP, 10._DP]
    point = [10._DP, 10._DP, 10._DP]
    
    write(output_unit, *) point_is_inside(box_origin, box_size, point)

end program test_point_is_inside
