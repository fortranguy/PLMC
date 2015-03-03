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
            real(DP), intent(inout) :: position(num_dimensions)
        end subroutine Abstract_Boxed_Spheres_fold

    end interface

    type, extends(Abstract_Boxed_Spheres), public :: Bulk_Boxed_Spheres
    contains
        procedure :: accessible_volume => Bulk_Boxed_Spheres_accessible_volume
        procedure :: is_inside => Bulk_Boxed_Spheres_is_inside
        procedure :: fold => Bulk_Boxed_Spheres_fold
    end type Bulk_Boxed_Spheres

    type, extends(Abstract_Boxed_Spheres), public :: Slab_Boxed_Spheres
    contains
        procedure :: accessible_volume => Slab_Boxed_Spheres_accessible_volume
        procedure :: is_inside => Slab_Boxed_Spheres_is_inside
        procedure :: fold => Slab_Boxed_Spheres_fold
    end type Slab_Boxed_Spheres

contains

!implementation Bulk_Boxed_Spheres

    pure function Bulk_Boxed_Spheres_accessible_volume(this) result(accessible_volume)
        class(Bulk_Boxed_Spheres), intent(in) :: this
        real(DP) :: accessible_volume

        accessible_volume = product(this%box%get_size())
    end function Bulk_Boxed_Spheres_accessible_volume

    pure function Bulk_Boxed_Spheres_is_inside(this, position) result(is_inside)
        class(Bulk_Boxed_Spheres), intent(in) :: this
        real(DP), intent(in) :: position(num_dimensions)
        logical :: is_inside

        is_inside = .true.
    end function Bulk_Boxed_Spheres_is_inside

    pure subroutine Bulk_Boxed_Spheres_fold(this, position)
        class(Bulk_Boxed_Spheres), intent(in) :: this
        real(DP), intent(inout) :: position(num_dimensions)

        position = modulo(position, this%box%get_size())
    end subroutine Bulk_Boxed_Spheres_fold

!end implementation Bulk_Boxed_Spheres

!implementation Slab_Boxed_Spheres

    pure function Slab_Boxed_Spheres_accessible_volume(this) result(accessible_volume)
        class(Slab_Boxed_Spheres), intent(in) :: this
        real(DP) :: accessible_volume

        accessible_volume = product(reshape(this%box%get_size(), [2])) * &
                            (this%box%get_height() - this%dip_spheres%get_diameter())
    end function Slab_Boxed_Spheres_accessible_volume

    pure function Slab_Boxed_Spheres_is_inside(this, position) result(is_inside)
        class(Slab_Boxed_Spheres), intent(in) :: this
        real(DP), intent(in) :: position(num_dimensions)
        logical :: is_inside

        is_inside = this%dip_spheres%get_diameter() / 2._DP < position(3) .and. &
                    position(3) < this%box%get_height() - this%dip_spheres%get_diameter() / 2._DP
    end function Slab_Boxed_Spheres_is_inside
        
    pure subroutine Slab_Boxed_Spheres_fold(this, position)
        class(Slab_Boxed_Spheres), intent(in) :: this
        real(DP), intent(inout) :: position(num_dimensions)

        position(1:2) = modulo(position(1:2), reshape(this%box%get_size(), [2]))
    end subroutine Slab_Boxed_Spheres_fold

!end implementation Slab_Boxed_Spheres

end module class_boxed_spheres
