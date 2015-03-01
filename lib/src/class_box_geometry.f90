module class_box_geometry

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions

implicit none

private

    type, abstract, public :: Abstract_Box_Geometry
    private
        real(DP) :: size(num_dimensions)
        integer :: wave(num_dimensions)
    contains
        procedure, non_overridable :: get_size => Abstract_Box_Geometry_get_size
        procedure(Abstract_Box_Geometry_get_height), deferred :: get_height
        procedure, non_overridable :: get_wave => Abstract_Box_Geometry_get_wave
        procedure(Abstract_Box_Geometry_vector_PBC), deferred :: vector_PBC
        procedure, non_overridable :: distance_PBC => Abstract_Box_Geometry_distance_PBC
    end type Abstract_Box_Geometry

    abstract interface

        pure function Abstract_Box_Geometry_get_height(this) result(get_height)
        import :: Abstract_Box_Geometry, DP
            class(Abstract_Box_Geometry), intent(in) :: this
            real(DP) :: get_height
        end function Abstract_Box_Geometry_get_height

        !> PBC: Periodic Boundary Conditions
        pure function Abstract_Box_Geometry_vector_PBC(this, position1, position2) result(vector_PBC)
        import :: Abstract_Box_Geometry, DP, num_dimensions
            class(Abstract_Box_Geometry), intent(in) :: this
            real(DP), intent(in) :: position1(:), position2(:)
            real(DP) :: vector_PBC(num_dimensions)
        end function Abstract_Box_Geometry_vector_PBC
        
    end interface

    type, extends(Abstract_Box_Geometry), public :: Bulk_Geometry
    contains
        procedure :: get_height => Bulk_Geometry_get_height
        procedure :: vector_PBC => Bulk_Geometry_vector_PBC
    end type Bulk_Geometry
    
contains

!implementation Abstract_Box_Geometry

    pure function Abstract_Box_Geometry_get_size(this) result(get_size)
        class(Abstract_Box_Geometry), intent(in) :: this
        real(DP) :: get_size(num_dimensions)

        get_size = this%size
    end function Abstract_Box_Geometry_get_size

    pure function Abstract_Box_Geometry_get_wave(this) result(get_wave)
        class(Abstract_Box_Geometry), intent(in) :: this
        integer :: get_wave(num_dimensions)

        get_wave = this%wave
    end function Abstract_Box_Geometry_get_wave

    pure function Abstract_Box_Geometry_distance_PBC(this, position1, position2) result(distance_PBC)
        class(Abstract_Box_Geometry), intent(in) :: this
        real(DP), intent(in) :: position1(:), position2(:)
        real(DP) :: distance_PBC

        distance_PBC = norm2(this%vector_PBC(position1, position2))
    end function Abstract_Box_Geometry_distance_PBC
    
!end implementation Abstract_Box_Geometry

!implementation Bulk_Geometry

    pure function Bulk_Geometry_get_height(this) result(get_height)
        class(Bulk_Geometry), intent(in) :: this
        real(DP) :: get_height

        get_height = this%size(3)
    end function Bulk_Geometry_get_height

    !> from SMAC, algorithm 2.5 & 2.6, p.91
    pure function Bulk_Geometry_vector_PBC(this, position1, position2) result(vector_PBC)
        class(Bulk_Geometry), intent(in) :: this
        real(DP), intent(in) :: position1(:), position2(:)
        real(DP) :: vector_PBC(num_dimensions)

        vector_PBC = modulo(position2 - position1, this%size)

        where(vector_PBC > this%size/2._DP)
            vector_PBC = vector_PBC - this%size
        end where
    end function Bulk_Geometry_vector_PBC

!end implementation Bulk_Geometry

end module class_box_geometry
