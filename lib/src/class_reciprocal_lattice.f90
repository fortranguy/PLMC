module class_reciprocal_lattice

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use module_error, only: warning_continue
use class_periodic_box, only: Abstract_Periodic_Box

implicit none

private

    type, abstract, public :: Abstract_Reciprocal_Lattice
    contains
        procedure(Abstract_Reciprocal_Lattice_set_reci_num), deferred :: set_reci_num
        procedure(Abstract_Reciprocal_Lattice_get_reci_num), deferred :: get_reci_num
    end type Abstract_Reciprocal_Lattice

    abstract interface
    
        subroutine Abstract_Reciprocal_Lattice_set_reci_num(this, periodic_box, reci_num)
        import :: num_dimensions, Abstract_Reciprocal_Lattice, Abstract_Periodic_Box
            class(Abstract_Reciprocal_Lattice), intent(out) :: this
            class(Abstract_Periodic_Box), intent(in) :: periodic_box
            integer, intent(in) :: reci_num(num_dimensions)
        end subroutine Abstract_Reciprocal_Lattice_set_reci_num

        pure function Abstract_Reciprocal_Lattice_get_reci_num(this) result(get_reci_num)
        import :: num_dimensions, Abstract_Reciprocal_Lattice
            class(Abstract_Reciprocal_Lattice), intent(in) :: this
            integer :: get_reci_num(num_dimensions)
        end function Abstract_Reciprocal_Lattice_get_reci_num
        
    end interface

    type, extends(Abstract_Reciprocal_Lattice), public :: Null_Reciprocal_Lattice
    contains
        procedure :: set_reci_num => Null_Reciprocal_Lattice_set_reci_num
        procedure :: get_reci_num => Null_Reciprocal_Lattice_get_reci_num
    end type Null_Reciprocal_Lattice

    type, extends(Abstract_Reciprocal_Lattice), public :: Reciprocal_Lattice
        private
        integer :: reci_num(num_dimensions)
    contains
        procedure :: set_reci_num => Reciprocal_Lattice_set_reci_num
        procedure, nopass, private :: check => Reciprocal_Lattice_check
        procedure :: get_reci_num => Reciprocal_Lattice_get_reci_num
    end type Reciprocal_Lattice
    
contains

!implementation Null_Reciprocal_Lattice

    subroutine Null_Reciprocal_Lattice_set_reci_num(this, periodic_box, reci_num)
        class(Null_Reciprocal_Lattice), intent(out) :: this
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        integer, intent(in) :: reci_num(num_dimensions)

    end subroutine Null_Reciprocal_Lattice_set_reci_num

    pure function Null_Reciprocal_Lattice_get_reci_num(this) result(get_reci_num)
        class(Null_Reciprocal_Lattice), intent(in) :: this
        integer :: get_reci_num(num_dimensions)

        get_reci_num = 0    
    end function Null_Reciprocal_Lattice_get_reci_num

!end implementation Reciprocal_Lattice

!implementation Null_Reciprocal_Lattice

    subroutine Reciprocal_Lattice_set_reci_num(this, periodic_box, reci_num)
        class(Reciprocal_Lattice), intent(out) :: this
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        integer, intent(in) :: reci_num(num_dimensions)

        call this%check(periodic_box%get_real_size(), reci_num)
        this%reci_num = reci_num
    end subroutine Reciprocal_Lattice_set_reci_num

    subroutine Reciprocal_Lattice_check(real_size, reci_num)
        real(DP), intent(in) :: real_size(num_dimensions)
        integer, intent(in) :: reci_num(num_dimensions)

        real(DP) :: real_zx_ratio, reci_zx_ratio
        
        if (reci_num(1) /= reci_num(2)) then
            call warning_continue("reci_num(1) and reci_num(2) are not equal.")
        end if

        real_zx_ratio = real_size(3) / real_size(1)
        reci_zx_ratio = real(reci_num(3), DP) / real(reci_num(1), DP)
        if (reci_zx_ratio < real_zx_ratio) then
            call warning_continue("reci z/x ratio is lower than real z/x ratio.")
        end if
    end subroutine Reciprocal_Lattice_check

    pure function Reciprocal_Lattice_get_reci_num(this) result(get_reci_num)
        class(Reciprocal_Lattice), intent(in) :: this
        integer :: get_reci_num(num_dimensions)

        get_reci_num = this%reci_num
    end function Reciprocal_Lattice_get_reci_num
    
!end implementation Reciprocal_Lattice

end module class_reciprocal_lattice
