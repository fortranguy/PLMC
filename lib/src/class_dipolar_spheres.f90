module class_dipolar_spheres

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_coordinates, only: increase_coordinates_size

implicit none

private

    type, abstract, public :: Abstract_Dipolar_Spheres
    private
        character(len=:), allocatable :: name
        real(DP) :: diameter
        integer ::  num
        real(DP), allocatable :: positions(:, :)
    contains
        procedure, non_overridable :: get_name => Abstract_Dipolar_Spheres_get_name
        procedure, non_overridable :: get_diameter => Abstract_Dipolar_Spheres_get_diameter
        procedure, non_overridable :: get_num => Abstract_Dipolar_Spheres_get_num
        procedure, non_overridable :: get_position => Abstract_Dipolar_Spheres_get_position
        procedure, non_overridable :: set_position => Abstract_Dipolar_Spheres_set_position
        procedure(Abstract_Dipolar_Spheres_get_moment_norm), deferred :: get_moment_norm
        procedure(Abstract_Dipolar_Spheres_get_moment), deferred :: get_moment
        procedure(Abstract_Dipolar_Spheres_set_moment), deferred :: set_moment
        procedure(Abstract_Dipolar_Spheres_remove), deferred :: remove
        procedure(Abstract_Dipolar_Spheres_add), deferred :: add
    end type Abstract_Dipolar_Spheres
    
    abstract interface
    
        pure function Abstract_Dipolar_Spheres_get_moment_norm(this) result(get_moment_norm)
        import :: Abstract_Dipolar_Spheres, DP
            class(Abstract_Dipolar_Spheres), intent(in) :: this
            real(DP) :: get_moment_norm
        end function Abstract_Dipolar_Spheres_get_moment_norm
        
        pure function Abstract_Dipolar_Spheres_get_moment(this, i_sphere) result(get_moment)
        import :: DP, num_dimensions, Abstract_Dipolar_Spheres
            class(Abstract_Dipolar_Spheres), intent(in) :: this
            integer, intent(in) :: i_sphere
            real(DP) :: get_moment(num_dimensions)
        end function Abstract_Dipolar_Spheres_get_moment
        
        pure subroutine Abstract_Dipolar_Spheres_set_moment(this, i_sphere, moment)
        import :: DP, num_dimensions, Abstract_Dipolar_Spheres
            class(Abstract_Dipolar_Spheres), intent(inout) :: this
            integer, intent(in) :: i_sphere
            real(DP), intent(in) :: moment(num_dimensions)
        end subroutine Abstract_Dipolar_Spheres_set_moment
        
        pure subroutine Abstract_Dipolar_Spheres_remove(this, i_sphere)
        import :: Abstract_Dipolar_Spheres
            class(Abstract_Dipolar_Spheres), intent(inout) :: this
            integer, intent(in) :: i_sphere            
        end subroutine Abstract_Dipolar_Spheres_remove
        
        pure subroutine Abstract_Dipolar_Spheres_add(this, position, moment)
        import :: DP, num_dimensions, Abstract_Dipolar_Spheres
            class(Abstract_Dipolar_Spheres), intent(inout) :: this
            real(DP), intent(in) :: position(num_dimensions)
            real(DP), intent(in) :: moment(num_dimensions)
        end subroutine Abstract_Dipolar_Spheres_add
    
    end interface
    
    type, extends(Abstract_Dipolar_Spheres), public :: Apolar_Spheres
    contains
        procedure :: get_moment_norm => Apolar_Spheres_get_moment_norm
        procedure :: get_moment => Apolar_Spheres_get_moment
        procedure :: set_moment => Apolar_Spheres_set_moment
        procedure :: remove => Apolar_Spheres_remove
        procedure :: add => Apolar_Spheres_add
    end type Apolar_Spheres
    
    type, extends(Abstract_Dipolar_Spheres), public :: Dipolar_Spheres
        real(DP) :: moment_norm
        real(DP), allocatable :: moments(:, :)
    contains
        procedure :: get_moment_norm => Dipolar_Spheres_get_moment_norm
        procedure :: get_moment => Dipolar_Spheres_get_moment
        procedure :: set_moment => Dipolar_Spheres_set_moment
        procedure :: remove => Dipolar_Spheres_remove
        procedure :: add => Dipolar_Spheres_add
    end type Dipolar_Spheres
    
contains

!implementation Abstract_Dipolar_Spheres

    pure function Abstract_Dipolar_Spheres_get_name(this) result(get_name)
        class(Abstract_Dipolar_Spheres), intent(in) :: this
        character(len=:), allocatable :: get_name
        
        get_name = this%name
    end function Abstract_Dipolar_Spheres_get_name
    
    pure function Abstract_Dipolar_Spheres_get_diameter(this) result(get_diameter)
        class(Abstract_Dipolar_Spheres), intent(in) :: this
        real(DP) :: get_diameter
        
        get_diameter = this%diameter
    end function Abstract_Dipolar_Spheres_get_diameter
    
    pure function Abstract_Dipolar_Spheres_get_num(this) result(get_num)
        class(Abstract_Dipolar_Spheres), intent(in) :: this
        integer :: get_num
        
        get_num = this%num
    end function Abstract_Dipolar_Spheres_get_num
    
    pure function Abstract_Dipolar_Spheres_get_position(this, i_sphere) result(get_position)
        class(Abstract_Dipolar_Spheres), intent(in) :: this
        integer, intent(in) :: i_sphere
        real(DP) :: get_position(num_dimensions)
        
        get_position = this%positions(:, i_sphere)
    end function Abstract_Dipolar_Spheres_get_position
    
    pure subroutine Abstract_Dipolar_Spheres_set_position(this, i_sphere, position)
        class(Abstract_Dipolar_Spheres), intent(inout) :: this
        integer, intent(in) :: i_sphere
        real(DP), intent(in) :: position(num_dimensions)
        
        this%positions(:, i_sphere) = position
    end subroutine Abstract_Dipolar_Spheres_set_position
    
!implementation Abstract_Dipolar_Spheres

!implementation Apolar_Spheres

    pure function Apolar_Spheres_get_moment_norm(this) result(get_moment_norm)
        class(Apolar_Spheres), intent(in) :: this
        real(DP) :: get_moment_norm
        
        get_moment_norm = 0._DP
    end function Apolar_Spheres_get_moment_norm
    
    pure function Apolar_Spheres_get_moment(this, i_sphere) result(get_moment)
        class(Apolar_Spheres), intent(in) :: this
        integer, intent(in) :: i_sphere
        real(DP) :: get_moment(num_dimensions)
        
        get_moment = 0._DP
    end function Apolar_Spheres_get_moment
    
    pure subroutine Apolar_Spheres_set_moment(this, i_sphere, moment)
        class(Apolar_Spheres), intent(inout) :: this
        integer, intent(in) :: i_sphere
        real(DP), intent(in) :: moment(num_dimensions)
        
    end subroutine Apolar_Spheres_set_moment
    
    pure subroutine Apolar_Spheres_remove(this, i_sphere)
        class(Apolar_Spheres), intent(inout) :: this
        integer, intent(in) :: i_sphere
        
        if (i_sphere < this%num) then
            call this%set_position(i_sphere, this%get_position(this%num))
        end if
        
        this%num = this%num - 1        
    end subroutine Apolar_Spheres_remove
    
    pure subroutine Apolar_Spheres_add(this, position, moment)
        class(Apolar_Spheres), intent(inout) :: this
        real(DP), intent(in) :: position(num_dimensions)
        real(DP), intent(in) :: moment(num_dimensions)
        
        this%num = this%num + 1
        
        if (size(this%positions) < this%num) then
            call increase_coordinates_size(this%positions)
        end if
        call this%set_position(this%num, position)
    end subroutine Apolar_Spheres_add
    
!end implementation Apolar_Spheres

!implementation Dipolar_Spheres

    pure function Dipolar_Spheres_get_moment_norm(this) result(get_moment_norm)
        class(Dipolar_Spheres), intent(in) :: this
        real(DP) :: get_moment_norm
        
        get_moment_norm = this%moment_norm
    end function Dipolar_Spheres_get_moment_norm
    
    pure function Dipolar_Spheres_get_moment(this, i_sphere) result(get_moment)
        class(Dipolar_Spheres), intent(in) :: this
        integer, intent(in) :: i_sphere
        real(DP) :: get_moment(num_dimensions)
        
        get_moment = this%moments(:, i_sphere)
    end function Dipolar_Spheres_get_moment
    
    pure subroutine Dipolar_Spheres_set_moment(this, i_sphere, moment)
        class(Dipolar_Spheres), intent(inout) :: this
        integer, intent(in) :: i_sphere
        real(DP), intent(in) :: moment(num_dimensions)
        
        this%moments(:, i_sphere) = moment
    end subroutine Dipolar_Spheres_set_moment
    
    pure subroutine Dipolar_Spheres_remove(this, i_sphere)
        class(Dipolar_Spheres), intent(inout) :: this
        integer, intent(in) :: i_sphere
        
        if (i_sphere < this%num) then
            call this%set_position(i_sphere, this%get_position(this%num))
            call this%set_moment(i_sphere, this%get_moment(this%num))
        end if
        
        this%num = this%num - 1        
    end subroutine Dipolar_Spheres_remove
    
    pure subroutine Dipolar_Spheres_add(this, position, moment)
        class(Dipolar_Spheres), intent(inout) :: this
        real(DP), intent(in) :: position(num_dimensions)
        real(DP), intent(in) :: moment(num_dimensions)
        
        this%num = this%num + 1
        
        if (size(this%positions) < this%num) then
            call increase_coordinates_size(this%positions)
        end if
        call this%set_position(this%num, position)
        
        if (size(this%moments) < this%num) then
            call increase_coordinates_size(this%moments)
        end if
        call this%set_moment(this%num, moment)
    end subroutine Dipolar_Spheres_add
    
!end implementation Dipolar_Spheres

end module class_dipolar_spheres
