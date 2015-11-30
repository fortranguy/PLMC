module class_ewald_reci_weight

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use class_periodic_box, only: Abstract_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use class_permittivity, only: Abstract_Permittivity
use class_ewald_convergence_parameter, only: Abstract_Ewald_Convergence_Parameter

private

    type, abstract, public :: Abstract_Ewald_Reci_Weight
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        integer :: reci_numbers(num_dimensions)
        real(DP) :: permittivity
        class(Abstract_Ewald_Convergence_Parameter), pointer :: alpha => null()
        real(DP), dimension(:, :, :), allocatable :: weight
    contains
        procedure :: construct => Abstract_Ewald_Reci_Weight_construct
        procedure :: destroy => Abstract_Ewald_Reci_Weight_destroy
        procedure :: reset => Abstract_Ewald_Reci_Weight_set
        procedure :: get => Abstract_Ewald_Reci_Weight_get
        procedure, private :: set => Abstract_Ewald_Reci_Weight_set
    end type Abstract_Ewald_Reci_Weight

    type, extends(Abstract_Ewald_Reci_Weight), public :: Concrete_Ewald_Reci_Weight

    end type Concrete_Ewald_Reci_Weight

    type, extends(Abstract_Ewald_Reci_Weight), public :: Null_Ewald_Reci_Weight
    contains
        procedure :: construct => Null_Ewald_Reci_Weight_construct
        procedure :: destroy => Null_Ewald_Reci_Weight_destroy
        procedure :: get => Null_Ewald_Reci_Weight_get
        procedure, private :: set => Null_Ewald_Reci_Weight_set
    end type Null_Ewald_Reci_Weight

contains

!implementation Abstract_Ewald_Reci_Weight

    subroutine Abstract_Ewald_Reci_Weight_construct(this, periodic_box, reciprocal_lattice, &
        permittivity, alpha)
        class(Abstract_Ewald_Reci_Weight), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Ewald_Convergence_Parameter), target, intent(in) :: alpha

        this%periodic_box => periodic_box
        this%reci_numbers = reciprocal_lattice%get_numbers()
        this%permittivity = permittivity%get()
        this%alpha => alpha
        allocate(this%weight(-this%reci_numbers(1):this%reci_numbers(1), &
                             -this%reci_numbers(2):this%reci_numbers(2), &
                             -this%reci_numbers(3):this%reci_numbers(3)))
        call this%set()
    end subroutine Abstract_Ewald_Reci_Weight_construct

    subroutine Abstract_Ewald_Reci_Weight_destroy(this)
        class(Abstract_Ewald_Reci_Weight), intent(inout) :: this

        if (allocated(this%weight)) deallocate(this%weight)
        this%alpha => null()
        this%periodic_box => null()
    end subroutine Abstract_Ewald_Reci_Weight_destroy

    pure subroutine Abstract_Ewald_Reci_Weight_set(this)
        class(Abstract_Ewald_Reci_Weight), intent(inout) :: this

        integer :: n_1, n_2, n_3
        real(DP), dimension(num_dimensions) :: box_size, wave_vector
        real(DP) :: alpha, wave_squared

        box_size = this%periodic_box%get_size()
        alpha = this%alpha%get()
        ! How about using symmetry to reduce memory?
        do n_3 = -this%reci_numbers(3), this%reci_numbers(3)
            wave_vector(3) = 2._DP*PI * real(n_3, DP) / box_size(3)
        do n_2 = -this%reci_numbers(2), this%reci_numbers(2)
            wave_vector(2) = 2._DP*PI * real(n_2, DP) / box_size(2)
        do n_1 = -this%reci_numbers(1), this%reci_numbers(1)
            wave_vector(1) = 2._DP*PI * real(n_1, DP) / box_size(1)
            if (n_1 /= 0 .or. n_2 /= 0 .or. n_3 /= 0) then
                wave_squared = dot_product(wave_vector, wave_vector)
                this%weight(n_1, n_2, n_3) = exp(-wave_squared/alpha**2/4._DP) / &
                    this%permittivity / product(box_size) / wave_squared
            else
                this%weight(n_1, n_2, n_3) = 0._DP
            end if
        end do
        end do
        end do
    end subroutine Abstract_Ewald_Reci_Weight_set

    !> \[
    !>      w_\alpha(\vec{k}) = \frac{e^{-k^2/4\alpha^2}}{\epsilon V k^2}
    !> \]
    pure real(DP) function Abstract_Ewald_Reci_Weight_get(this, n_1, n_2, n_3) result(weight)
        class(Abstract_Ewald_Reci_Weight), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3

        weight = this%weight(n_1, n_2, n_3)
    end function Abstract_Ewald_Reci_Weight_get

!end implementation Abstract_Ewald_Reci_Weight

!implementation Null_Ewald_Reci_Weight

    subroutine Null_Ewald_Reci_Weight_construct(this, periodic_box, reciprocal_lattice, &
        permittivity, alpha)
        class(Null_Ewald_Reci_Weight), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Ewald_Convergence_Parameter), target, intent(in) :: alpha
    end subroutine Null_Ewald_Reci_Weight_construct

    subroutine Null_Ewald_Reci_Weight_destroy(this)
        class(Null_Ewald_Reci_Weight), intent(inout) :: this
    end subroutine Null_Ewald_Reci_Weight_destroy

    pure subroutine Null_Ewald_Reci_Weight_set(this)
        class(Null_Ewald_Reci_Weight), intent(inout) :: this
    end subroutine Null_Ewald_Reci_Weight_set

    pure real(DP) function Null_Ewald_Reci_Weight_get(this, n_1, n_2, n_3) result(weight)
        class(Null_Ewald_Reci_Weight), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3
        weight = 0._DP
    end function Null_Ewald_Reci_Weight_get

!end implementation Null_Ewald_Reci_Weight

end module class_ewald_reci_weight
