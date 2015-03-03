module class_boxed_spheres

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_box_geometry, only: Abstract_Box_Geometry
use class_dipolar_spheres, only: Abstract_Dipolar_Spheres

implicit none

private

    type, abstract, public :: Abstract_Boxed_Spheres
    private
        class(Abstract_Box_Geometry), pointer :: box => null()
        class(Abstract_Dipolar_Spheres), pointer :: dip_spheres => null()
    contains
        procedure(Abstract_Boxed_Spheres_accessible_volume), deferred :: accessible_volume
        procedure(Abstract_Boxed_Spheres_is_inside), deferred :: is_inside
        procedure(Abstract_Boxed_Spheres_fold), deferred :: fold
    end type Abstract_Boxed_Spheres

    abstract interface

        pure function Abstract_Boxed_Spheres_accessible_volume(this) result(accessible_volume)
        import :: DP, Abstract_Boxed_Spheres
            class(Abstract_Boxed_Spheres), intent(in) :: this
            real(DP) :: accessible_volume
        end function Abstract_Boxed_Spheres_accessible_volume

        pure function Abstract_Boxed_Spheres_is_inside(this, position) result(is_inside)
        import :: DP, num_dimensions, Abstract_Boxed_Spheres
            class(Abstract_Boxed_Spheres), intent(in) :: this
            real(DP), intent(in) :: position(num_dimensions)
            logical :: is_inside
        end function Abstract_Boxed_Spheres_is_inside

        pure subroutine Abstract_Boxed_Spheres_fold(this, position)
        import :: DP, num_dimensions, Abstract_Boxed_Spheres
            class(Abstract_Boxed_Spheres), intent(in) :: this
            real(DP), intent(in) :: position(num_dimensions)
        end subroutine Abstract_Boxed_Spheres_fold

    end interface

end module class_boxed_spheres
