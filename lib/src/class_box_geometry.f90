module class_box_geometry

use, intrinsic :: iso_fortran_env, only: DP => REAL64, error_unit
use data_geometry, only: num_dimensions
use json_module, only: json_file
use module_data, only: test_data_found

implicit none

private

    type, abstract, public :: Abstract_Box_Geometry
    private
        real(DP) :: size(num_dimensions)
        integer :: wave(num_dimensions)
    contains
        procedure(Abstract_Box_Geometry_set), deferred :: set
        procedure, private :: set_size => Abstract_Box_Geometry_set_size
        procedure, private, nopass :: check_size => Abstract_Box_Geometry_check_size
        procedure, private :: set_wave => Abstract_Box_Geometry_set_wave
        procedure, private, nopass :: check_wave => Abstract_Box_Geometry_check_wave
        
        procedure, non_overridable :: get_size => Abstract_Box_Geometry_get_size
        procedure(Abstract_Box_Geometry_get_height), deferred :: get_height
        procedure, non_overridable :: get_wave => Abstract_Box_Geometry_get_wave
        
        procedure(Abstract_Box_Geometry_vector_PBC), deferred :: vector_PBC
        procedure, non_overridable :: distance_PBC => Abstract_Box_Geometry_distance_PBC
    end type Abstract_Box_Geometry

    abstract interface

        subroutine Abstract_Box_Geometry_set(this, input_data)
        import :: Abstract_Box_Geometry, json_file
            class(Abstract_Box_Geometry), intent(out) :: this
            type(json_file), intent(inout) :: input_data
        end subroutine

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
        procedure :: set => Bulk_Geometry_set
        procedure :: get_height => Bulk_Geometry_get_height
        procedure :: vector_PBC => Bulk_Geometry_vector_PBC
    end type Bulk_Geometry

    type, extends(Abstract_Box_Geometry), public :: Slab_Geometry
    private
        real(DP) :: height
    contains
        procedure :: set => Slab_Geometry_set
        procedure, private :: set_height => Slab_Geometry_set_height
        procedure :: get_height => Slab_Geometry_get_height
        procedure :: vector_PBC => Slab_Geometry_vector_PBC
    end type Slab_Geometry
    
contains

!implementation Abstract_Box_Geometry

    subroutine Abstract_Box_Geometry_set_size(this, input_data, size_dimension)
        class(Abstract_Box_Geometry), intent(inout) :: this
        type(json_file), intent(inout) :: input_data
        integer, intent(in) :: size_dimension

        character(len=:), allocatable :: data_field
        logical :: found
        real(DP), allocatable :: Box_size(:) ! workaround: this%size can't be set.

        data_field = "Box Geometry.size"
        call input_data%get(data_field, Box_size, found)
        call test_data_found(data_field, found)

        call this%check_size(data_field, Box_size, size_dimension)

        this%size(1:size_dimension) = Box_size
    end subroutine Abstract_Box_Geometry_set_size

    subroutine Abstract_Box_Geometry_check_size(size_field, size_value, size_dimension)
        character(len=*), intent(in) :: size_field
        real(DP), intent(in) :: size_value(:)
        integer, intent(in) :: size_dimension

        if (size(size_value) /= size_dimension) then
            write(error_unit, *) size_field, " dimension is not ", size_dimension, "."
            error stop
        end if
    end subroutine Abstract_Box_Geometry_check_size

    subroutine Abstract_Box_Geometry_set_wave(this, input_data, wave_dimension)
        class(Abstract_Box_Geometry), intent(inout) :: this
        type(json_file), intent(inout) :: input_data
        integer, intent(in) :: wave_dimension

        character(len=:), allocatable :: data_field
        logical :: found
        integer, allocatable :: Box_wave(:) ! workaround: this%wave can't be set.

        data_field = "Box Geometry.wave"
        call input_data%get(data_field, Box_wave, found)
        call test_data_found(data_field, found)

        call this%check_wave(data_field, Box_wave, wave_dimension)

        this%wave(1:wave_dimension) = Box_wave
    end subroutine Abstract_Box_Geometry_set_wave

    subroutine Abstract_Box_Geometry_check_wave(wave_field, wave_value, wave_dimension)
        character(len=*), intent(in) :: wave_field
        integer, intent(in) :: wave_value(:)
        integer, intent(in) :: wave_dimension

        if (size(wave_value) /= wave_dimension) then
            write(error_unit, *) wave_field, " dimension is not ", wave_dimension, "."
            error stop
        end if
    end subroutine Abstract_Box_Geometry_check_wave

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

    subroutine Bulk_Geometry_set(this, input_data)
        class(Bulk_Geometry), intent(out) :: this
        type(json_file), intent(inout) :: input_data

        call this%set_size(input_data, num_dimensions)
        call this%set_wave(input_data, num_dimensions)        
    end subroutine Bulk_Geometry_set

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

!implementation Slab_Geometry

    subroutine Slab_Geometry_set(this, input_data)
        class(Slab_Geometry), intent(out) :: this
        type(json_file), intent(inout) :: input_data

        call this%set_size(input_data, 2)
        call this%set_wave(input_data, 2)
        call this%set_height(input_data)        
    end subroutine Slab_Geometry_set

    subroutine Slab_Geometry_set_height(this, input_data)
        class(Slab_Geometry), intent(inout) :: this
        type(json_file), intent(inout) :: input_data

        character(len=:), allocatable :: data_field
        logical :: found

        data_field = "Box Geometry.height"
        call input_data%get(data_field, this%height, found)
        call test_data_found(data_field, found)        
    end subroutine Slab_Geometry_set_height

    pure function Slab_Geometry_get_height(this) result(get_height)
        class(Slab_Geometry), intent(in) :: this
        real(DP) :: get_height

        get_height = this%height
    end function Slab_Geometry_get_height

    pure function Slab_Geometry_vector_PBC(this, position1, position2) result(vector_PBC)
        class(Slab_Geometry), intent(in) :: this
        real(DP), intent(in) :: position1(:), position2(:)
        real(DP) :: vector_PBC(num_dimensions)

        vector_PBC(1:2) = modulo(position2(1:2) - position1(1:2), this%size(1:2))

        where(vector_PBC(1:2) > this%size(1:2)/2._DP)
            vector_PBC(1:2) = vector_PBC(1:2) - this%size(1:2)
        end where

        vector_PBC(3) = position2(3) - position1(3)
    end function Slab_Geometry_vector_PBC

!end implementation Slab_Geometry

end module class_box_geometry
