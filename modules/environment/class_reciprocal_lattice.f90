module class_reciprocal_lattice

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_errors, only: warning_continue, error_exit
use procedures_checks, only: check_3d_array, check_positive
use class_periodic_box, only: Abstract_Periodic_Box

implicit none

private

    type, abstract, public :: Abstract_Reciprocal_Lattice
    private
        integer :: numbers(num_dimensions)
        class(Abstract_Periodic_Box), pointer :: periodic_box
    contains
        procedure :: construct => Abstract_Reciprocal_Lattice_construct
        procedure, private :: check => Abstract_Reciprocal_Lattice_check
        procedure :: destroy => Abstract_Reciprocal_Lattice_destroy
        procedure :: get_numbers => Abstract_Reciprocal_Lattice_get_numbers
    end type Abstract_Reciprocal_Lattice

    type, extends(Abstract_Reciprocal_Lattice), public :: Concrete_Reciprocal_Lattice

    end type Concrete_Reciprocal_Lattice

    type, extends(Abstract_Reciprocal_Lattice), public :: Null_Reciprocal_Lattice
    contains
        procedure :: construct => Null_Reciprocal_Lattice_construct
        procedure :: destroy => Null_Reciprocal_Lattice_destroy
        procedure :: get_numbers => Null_Reciprocal_Lattice_get_numbers
    end type Null_Reciprocal_Lattice

contains

!implementation Abstract_Reciprocal_Lattice

    subroutine Abstract_Reciprocal_Lattice_construct(this, periodic_box, numbers)
        class(Abstract_Reciprocal_Lattice), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        integer, intent(in) :: numbers(:)

        this%periodic_box => periodic_box
        call this%check(numbers)
        this%numbers = numbers
    end subroutine Abstract_Reciprocal_Lattice_construct

    subroutine Abstract_Reciprocal_Lattice_check(this, numbers)
        class(Abstract_Reciprocal_Lattice), intent(in) :: this
        integer, intent(in) :: numbers(:)

        real(DP) :: real_size(num_dimensions)
        real(DP) :: real_zx_ratio, reci_zx_ratio

        call check_3d_array("Abstract_Reciprocal_Lattice", "numbers", numbers)
        call check_positive("Abstract_Reciprocal_Lattice", "numbers", numbers)
        if (numbers(1) /= numbers(2)) then
            call warning_continue("Abstract_Reciprocal_Lattice: "//&
                "numbers(1) and numbers(2) are not equal.")
        end if

        real_size = this%periodic_box%get_size()
        real_zx_ratio = real_size(3) / real_size(1)
        reci_zx_ratio = real(numbers(3), DP) / real(numbers(1), DP)
        if (reci_zx_ratio < real_zx_ratio) then
            call warning_continue("Abstract_Reciprocal_Lattice: "//&
                "reci z/x ratio is lower than real z/x ratio.")
        end if
    end subroutine Abstract_Reciprocal_Lattice_check

    subroutine Abstract_Reciprocal_Lattice_destroy(this)
        class(Abstract_Reciprocal_Lattice), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Abstract_Reciprocal_Lattice_destroy

    pure function Abstract_Reciprocal_Lattice_get_numbers(this) result(numbers)
        class(Abstract_Reciprocal_Lattice), intent(in) :: this
        integer :: numbers(num_dimensions)

        numbers = this%numbers
    end function Abstract_Reciprocal_Lattice_get_numbers

!end implementation Abstract_Reciprocal_Lattice

!implementation Null_Reciprocal_Lattice

    subroutine Null_Reciprocal_Lattice_construct(this, periodic_box, numbers)
        class(Null_Reciprocal_Lattice), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        integer, intent(in) :: numbers(:)
    end subroutine Null_Reciprocal_Lattice_construct

    subroutine Null_Reciprocal_Lattice_destroy(this)
        class(Null_Reciprocal_Lattice), intent(inout) :: this
    end subroutine Null_Reciprocal_Lattice_destroy

    pure function Null_Reciprocal_Lattice_get_numbers(this) result(numbers)
        class(Null_Reciprocal_Lattice), intent(in) :: this
        integer :: numbers(num_dimensions)
        numbers = 0
    end function Null_Reciprocal_Lattice_get_numbers

!end implementation Null_Reciprocal_Lattice

end module class_reciprocal_lattice
