module class_reciprocal_lattice

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_errors, only: warning_continue, error_exit
use class_periodic_box, only: Abstract_Periodic_Box

implicit none

private

    type, abstract, public :: Abstract_Reciprocal_Lattice
    contains
        procedure(Abstract_Reciprocal_Lattice_construct), deferred :: construct
        procedure(Abstract_Reciprocal_Lattice_destroy), deferred :: destroy
        procedure(Abstract_Reciprocal_Lattice_get_reci_num), deferred :: get_reci_num
    end type Abstract_Reciprocal_Lattice

    abstract interface
    
        subroutine Abstract_Reciprocal_Lattice_construct(this, periodic_box, reci_num)
        import :: Abstract_Reciprocal_Lattice, Abstract_Periodic_Box
            class(Abstract_Reciprocal_Lattice), intent(out) :: this
            class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
            integer, intent(in) :: reci_num(:)
        end subroutine Abstract_Reciprocal_Lattice_construct
        
        subroutine Abstract_Reciprocal_Lattice_destroy(this)
        import :: Abstract_Reciprocal_Lattice
            class(Abstract_Reciprocal_Lattice), intent(inout) :: this
        end subroutine Abstract_Reciprocal_Lattice_destroy

        pure function Abstract_Reciprocal_Lattice_get_reci_num(this) result(reci_num)
        import :: num_dimensions, Abstract_Reciprocal_Lattice
            class(Abstract_Reciprocal_Lattice), intent(in) :: this
            integer :: reci_num(num_dimensions)
        end function Abstract_Reciprocal_Lattice_get_reci_num
        
    end interface

    type, extends(Abstract_Reciprocal_Lattice), public :: Null_Reciprocal_Lattice
    contains
        procedure :: construct => Null_Reciprocal_Lattice_construct
        procedure :: destroy => Null_Reciprocal_Lattice_destroy
        procedure :: get_reci_num => Null_Reciprocal_Lattice_get_reci_num
    end type Null_Reciprocal_Lattice

    type, extends(Abstract_Reciprocal_Lattice), public :: Concrete_Reciprocal_Lattice
        private
        class(Abstract_Periodic_Box), pointer :: periodic_box
        integer :: reci_num(num_dimensions)
    contains
        procedure :: construct => Concrete_Reciprocal_Lattice_construct
        procedure, private :: check => Concrete_Reciprocal_Lattice_check
        procedure :: destroy => Concrete_Reciprocal_Lattice_destroy
        procedure :: get_reci_num => Concrete_Reciprocal_Lattice_get_reci_num
    end type Concrete_Reciprocal_Lattice
    
contains

!implementation Null_Reciprocal_Lattice

    subroutine Null_Reciprocal_Lattice_construct(this, periodic_box, reci_num)
        class(Null_Reciprocal_Lattice), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        integer, intent(in) :: reci_num(:)

    end subroutine Null_Reciprocal_Lattice_construct
    
    subroutine Null_Reciprocal_Lattice_destroy(this)
        class(Null_Reciprocal_Lattice), intent(inout) :: this
        
    end subroutine Null_Reciprocal_Lattice_destroy

    pure function Null_Reciprocal_Lattice_get_reci_num(this) result(reci_num)
        class(Null_Reciprocal_Lattice), intent(in) :: this
        integer :: reci_num(num_dimensions)

        reci_num = 0
    end function Null_Reciprocal_Lattice_get_reci_num

!end implementation Null_Reciprocal_Lattice

!implementation Concrete_Reciprocal_Lattice

    subroutine Concrete_Reciprocal_Lattice_construct(this, periodic_box, reci_num)
        class(Concrete_Reciprocal_Lattice), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        integer, intent(in) :: reci_num(:)
        
        this%periodic_box => periodic_box
        call this%check(reci_num)
        this%reci_num = reci_num
    end subroutine Concrete_Reciprocal_Lattice_construct

    subroutine Concrete_Reciprocal_Lattice_check(this, reci_num)
        class(Concrete_Reciprocal_Lattice), intent(in) :: this
        integer, intent(in) :: reci_num(:)
        
        real(DP) :: real_size(num_dimensions)
        real(DP) :: real_zx_ratio, reci_zx_ratio
        
        if (size(reci_num) /= num_dimensions) then
            call error_exit("Concrete_Reciprocal_Lattice: wrong number of dimensions (size).")
        end if
        if (reci_num(1) /= reci_num(2)) then
            call warning_continue("Concrete_Reciprocal_Lattice: "//&
                "reci_num(1) and reci_num(2) are not equal.")
        end if
        
        real_size = this%periodic_box%get_real_size()
        real_zx_ratio = real_size(3) / real_size(1)
        reci_zx_ratio = real(reci_num(3), DP) / real(reci_num(1), DP)
        if (reci_zx_ratio < real_zx_ratio) then
            call warning_continue("Concrete_Reciprocal_Lattice: "//&
                "reci z/x ratio is lower than real z/x ratio.")
        end if
    end subroutine Concrete_Reciprocal_Lattice_check
    
    subroutine Concrete_Reciprocal_Lattice_destroy(this)
        class(Concrete_Reciprocal_Lattice), intent(inout) :: this
        
        this%periodic_box => null()
    end subroutine Concrete_Reciprocal_Lattice_destroy

    pure function Concrete_Reciprocal_Lattice_get_reci_num(this) result(reci_num)
        class(Concrete_Reciprocal_Lattice), intent(in) :: this
        integer :: reci_num(num_dimensions)

        reci_num = this%reci_num
    end function Concrete_Reciprocal_Lattice_get_reci_num
    
!end implementation Concrete_Reciprocal_Lattice

end module class_reciprocal_lattice
