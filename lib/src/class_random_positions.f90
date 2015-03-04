module class_random_positions

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_box_geometry, only: Abstract_Box_Geometry
use class_dipolar_spheres, only: Abstract_Dipolar_Spheres

implicit none

private

    type, abstract, public :: Abstract_Random_Positions
    private
        class(Abstract_Box_Geometry), pointer :: box => null()
        class(Abstract_Dipolar_Spheres), pointer :: dip_spheres => null()
        real(DP) :: delta
    contains
        procedure :: construct => Abstract_Random_Positions_construct
        procedure :: destroy => Abstract_Random_Positions_destroy
        
        procedure(Abstract_Random_Positions_position), deferred :: position
        procedure :: move => Abstract_Random_Positions_move
    end type Abstract_Random_Positions
    
    abstract interface
    
        function Abstract_Random_Positions_position(this) result(position)
        import :: DP, num_dimensions, Abstract_Random_Positions
            class(Abstract_Random_Positions), intent(in) :: this
            real(DP) :: position(num_dimensions)
        end function Abstract_Random_Positions_position
        
    end interface
    
    type, extends(Abstract_Random_Positions), public :: Bulk_Random_Positions
    contains
        procedure :: position => Bulk_Random_Positions_position
    end type Bulk_Random_Positions
    
    type, extends(Abstract_Random_Positions), public :: Slab_Random_Positions
    contains
        procedure :: position => Slab_Random_Positions_position
    end type Slab_Random_Positions
    
contains

!implementation Abstract_Random_Positions

    subroutine Abstract_Random_Positions_construct(this, box, dip_spheres, delta)
        class(Abstract_Random_Positions), intent(out) :: this
        class(Abstract_Box_Geometry), target, intent(in) :: box
        class(Abstract_Dipolar_Spheres), target, intent(in) :: dip_spheres
        real(DP), intent(in) :: delta
        
        this%box => box
        this%dip_spheres => dip_spheres
        this%delta = delta
    end subroutine Abstract_Random_Positions_construct
    
    subroutine Abstract_Random_Positions_destroy(this)
        class(Abstract_Random_Positions), intent(inout) :: this
        
        this%dip_spheres => null()
        this%box => null()
    end subroutine Abstract_Random_Positions_destroy

    function Abstract_Random_Positions_move(this, i_sphere) result(position)
        class(Abstract_Random_Positions), intent(in) :: this
        integer, intent(in) :: i_sphere
        real(DP) :: position(num_dimensions)
        
        real(DP) :: random_vector(num_dimensions)
        
        position = this%dip_spheres%get_position(i_sphere)
        call random_number(random_vector)
        position = position + (random_vector - 0.5_DP) * this%delta
    end function Abstract_Random_Positions_move
    
!end implementation Abstract_Random_Positions

!implementation Bulk_Random_Positions

    function Bulk_Random_Positions_position(this) result(position)
        class(Bulk_Random_Positions), intent(in) :: this
        real(DP) :: position(num_dimensions)
        
        real(DP) :: random_vector(num_dimensions)
        
        call random_number(random_vector)
        position = this%box%get_size() * random_vector
    end function Bulk_Random_Positions_position

!end implementation Bulk_Random_Positions

!implementation Slab_Random_Positions

    function Slab_Random_Positions_position(this) result(position)
        class(Slab_Random_Positions), intent(in) :: this
        real(DP) :: position(num_dimensions)
        
        real(DP) :: random_vector(num_dimensions)
        
        call random_number(random_vector)
        position(1:2) = reshape(this%box%get_size(), [2]) * random_vector(1:2)
        position(3) = (this%box%get_height() - this%dip_spheres%get_diameter()) * &
                      random_vector(3) + this%dip_spheres%get_diameter() / 2._DP
    end function Slab_Random_Positions_position

!end implementation Slab_Random_Positions

end module class_random_positions
