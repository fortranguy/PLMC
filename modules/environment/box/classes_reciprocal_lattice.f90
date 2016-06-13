module classes_reciprocal_lattice

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_errors, only: warning_continue, error_exit
use procedures_checks, only: check_array_size, check_positive
use classes_periodic_box, only: Abstract_Periodic_Box

implicit none

private

    type, abstract, public :: Abstract_Reciprocal_Lattice
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        integer :: numbers(num_dimensions) = 0
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: is_sparse_in_z => Abstract_is_sparse_in_z
        procedure :: get_numbers => Abstract_get_numbers
    end type Abstract_Reciprocal_Lattice

    type, extends(Abstract_Reciprocal_Lattice), public :: Concrete_Reciprocal_Lattice

    end type Concrete_Reciprocal_Lattice

    type, extends(Abstract_Reciprocal_Lattice), public :: Null_Reciprocal_Lattice
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: is_sparse_in_z => Null_is_sparse_in_z
        procedure :: get_numbers => Null_get_numbers
    end type Null_Reciprocal_Lattice

contains

!implementation Abstract_Reciprocal_Lattice

    subroutine Abstract_construct(this, periodic_box, numbers)
        class(Abstract_Reciprocal_Lattice), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        integer, intent(in) :: numbers(:)

        this%periodic_box => periodic_box
        call check_array_size("Abstract_Reciprocal_Lattice", "numbers", numbers, num_dimensions)
        call check_positive("Abstract_Reciprocal_Lattice", "numbers", numbers)
        if (numbers(1) /= numbers(2)) then
            call warning_continue("Abstract_Reciprocal_Lattice: "//&
                "numbers(1) and numbers(2) are not equal.")
        end if
        this%numbers = numbers
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Reciprocal_Lattice), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Abstract_destroy

    !> This function tells if the density of wave vectors is too sparse in \( z \) direction.
    pure logical function Abstract_is_sparse_in_z(this) result(is_sparse_in_z)
        class(Abstract_Reciprocal_Lattice), intent(in) :: this

        real(DP) :: real_size(num_dimensions)
        real(DP) :: real_zx_ratio, reci_zx_ratio

        real_size = this%periodic_box%get_size()
        real_zx_ratio = real_size(3) / real_size(1)
        reci_zx_ratio = real(this%numbers(3), DP) / real(this%numbers(1), DP)
        if (reci_zx_ratio < real_zx_ratio) then
            is_sparse_in_z = .true.
        else
            is_sparse_in_z = .false.
        end if
    end function Abstract_is_sparse_in_z

    pure function Abstract_get_numbers(this) result(numbers)
        class(Abstract_Reciprocal_Lattice), intent(in) :: this
        integer :: numbers(num_dimensions)

        numbers = this%numbers
    end function Abstract_get_numbers

!end implementation Abstract_Reciprocal_Lattice

!implementation Null_Reciprocal_Lattice

    subroutine Null_construct(this, periodic_box, numbers)
        class(Null_Reciprocal_Lattice), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        integer, intent(in) :: numbers(:)
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Reciprocal_Lattice), intent(inout) :: this
    end subroutine Null_destroy

    pure logical function Null_is_sparse_in_z(this) result(is_sparse_in_z)
        class(Null_Reciprocal_Lattice), intent(in) :: this
        is_sparse_in_z = .false.
    end function Null_is_sparse_in_z

    pure function Null_get_numbers(this) result(numbers)
        class(Null_Reciprocal_Lattice), intent(in) :: this
        integer :: numbers(num_dimensions)
        numbers = 0
    end function Null_get_numbers

!end implementation Null_Reciprocal_Lattice

end module classes_reciprocal_lattice
