module class_parallelepiped_domain

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_periodic_box, only: Abstract_Periodic_Box
use procedures_checks, only: check_3d_array, check_positive

implicit none

private

    type, abstract, public :: Abstract_Parallelepiped_Domain
    private
        real(DP) :: origin(num_dimensions), size(num_dimensions)
        class(Abstract_Periodic_Box), pointer :: periodic_box
    contains
        procedure :: construct => Abstract_Parallelepiped_Domain_construct
        !procedure :: destroy => Abstract_Parallelepiped_Domain_destroy
        !procedure :: get_volume => Abstract_Parallelepiped_Domain_get_volume
        !procedure :: is_inside => Abstract_Parallelepiped_Domain_is_inside
    end type Abstract_Parallelepiped_Domain
    
contains

!implementation Abstract_Parallelepiped_Domain

    subroutine Abstract_Parallelepiped_Domain_construct(this, periodic_box, origin, size)
        class(Abstract_Parallelepiped_Domain), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: origin(:), size(:)
        
        this%periodic_box => periodic_box
        call check_3d_array("Abstract_Parallelepiped_Domain", "origin", origin)
        this%origin = this%origin
        call check_3d_array("Abstract_Parallelepiped_Domain", "size", size)
        call check_positive("Abstract_Parallelepiped_Domain", "size", size)
        this%size = size
    end subroutine Abstract_Parallelepiped_Domain_construct
    
!end implementation Abstract_Parallelepiped_Domain

end module class_parallelepiped_domain
