module procedures_parallelepiped_domain_micro

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions

implicit none

    abstract interface

        pure logical function abstract_is_inside(box_size, point_from_corner)
        import :: DP
            real(DP), intent(in) :: box_size(:), point_from_corner(:)
        end function abstract_is_inside

    end interface

private
public :: point_is_inside_box, is_inside_included, is_inside_excluded

contains

    pure logical function point_is_inside_box(box_origin, box_size, point, is_inside)
        real(DP), intent(in) :: box_origin(:), box_size(:), point(:)
        procedure(abstract_is_inside) :: is_inside

        real(DP), dimension(num_dimensions) :: box_corner, point_from_corner

        box_corner = box_origin - box_size/2._DP
        point_from_corner = point - box_corner
        point_is_inside_box = is_inside(box_size, point_from_corner)
    end function point_is_inside_box

    pure logical function is_inside_included(box_size, point_from_corner) result(is_inside)
        real(DP), intent(in) :: box_size(:), point_from_corner(:)

        is_inside = all(0._DP <= point_from_corner .and. point_from_corner <= box_size)
    end function is_inside_included

    pure logical function is_inside_excluded(box_size, point_from_corner) result(is_inside)
        real(DP), intent(in) :: box_size(:), point_from_corner(:)

        is_inside = all(0._DP <= point_from_corner .and. point_from_corner < box_size)
    end function is_inside_excluded

end module procedures_parallelepiped_domain_micro
