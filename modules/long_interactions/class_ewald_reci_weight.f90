module class_ewald_reci_weight

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use class_periodic_box, only: Abstract_Periodic_Box
use class_permittivity, only: Abstract_Permittivity
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
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
        procedure :: get_reci_numbers => Abstract_Ewald_Reci_Weight_get_reci_numbers
        procedure :: get => Abstract_Ewald_Reci_Weight_get
        procedure, private :: set => Abstract_Ewald_Reci_Weight_set
    end type Abstract_Ewald_Reci_Weight

    type, extends(Abstract_Ewald_Reci_Weight), public :: Concrete_Ewald_Reci_Weight

    end type Concrete_Ewald_Reci_Weight

    type, extends(Abstract_Ewald_Reci_Weight), public :: Null_Ewald_Reci_Weight
    contains
        procedure :: construct => Null_Ewald_Reci_Weight_construct
        procedure :: destroy => Null_Ewald_Reci_Weight_destroy
        procedure :: reset => Null_Ewald_Reci_Weight_set
        procedure :: get_reci_numbers => Null_Ewald_Reci_Weight_get_reci_numbers
        procedure :: get => Null_Ewald_Reci_Weight_get
    end type Null_Ewald_Reci_Weight

contains

!implementation Abstract_Ewald_Reci_Weight

    subroutine Abstract_Ewald_Reci_Weight_construct(this, periodic_box, permittivity, &
        reciprocal_lattice, alpha)
        class(Abstract_Ewald_Reci_Weight), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Ewald_Convergence_Parameter), target, intent(in) :: alpha

        this%periodic_box => periodic_box
        this%reci_numbers = reciprocal_lattice%get_numbers()
        this%permittivity = permittivity%get()
        this%alpha => alpha
        allocate(this%weight(0:this%reci_numbers(1), 0:this%reci_numbers(2), &
                             0:this%reci_numbers(3)))
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
        do n_3 = 0, this%reci_numbers(3)
            wave_vector(3) = 2._DP*PI * real(n_3, DP) / box_size(3)
        do n_2 = 0, this%reci_numbers(2)
            wave_vector(2) = 2._DP*PI * real(n_2, DP) / box_size(2)
        do n_1 = 0, this%reci_numbers(1)
            wave_vector(1) = 2._DP*PI * real(n_1, DP) / box_size(1)
            if (n_1 == 0 .and. n_2 == 0 .and. n_3 == 0) then
                this%weight(n_1, n_2, n_3) = 0._DP
            else
                wave_squared = dot_product(wave_vector, wave_vector)
                this%weight(n_1, n_2, n_3) = exp(-wave_squared/alpha**2/4._DP) / &
                    this%permittivity / product(box_size) / wave_squared
            end if
        end do
        end do
        end do
    end subroutine Abstract_Ewald_Reci_Weight_set

    pure function Abstract_Ewald_Reci_Weight_get_reci_numbers(this) result(reci_numbers)
        class(Abstract_Ewald_Reci_Weight), intent(in) :: this
        integer :: reci_numbers(num_dimensions)

        reci_numbers = this%reci_numbers
    end function Abstract_Ewald_Reci_Weight_get_reci_numbers

    !> \[
    !>      w_\alpha(\vec{k}) =
    !>          \begin{cases}
    !>              0                                           & \text{if } \vec{k} = \vec{0} \\
    !>              \frac{e^{-k^2/4\alpha^2}}{\epsilon V k^2}   & \text{else}
    !>          \end{cases}
    !> \]
    pure real(DP) function Abstract_Ewald_Reci_Weight_get(this, n_1, n_2, n_3) result(weight)
        class(Abstract_Ewald_Reci_Weight), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3

        weight = this%weight(abs(n_1), abs(n_2), abs(n_3))
    end function Abstract_Ewald_Reci_Weight_get

!end implementation Abstract_Ewald_Reci_Weight

!implementation Null_Ewald_Reci_Weight

    subroutine Null_Ewald_Reci_Weight_construct(this, periodic_box, permittivity, &
        reciprocal_lattice, alpha)
        class(Null_Ewald_Reci_Weight), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Ewald_Convergence_Parameter), target, intent(in) :: alpha
    end subroutine Null_Ewald_Reci_Weight_construct

    subroutine Null_Ewald_Reci_Weight_destroy(this)
        class(Null_Ewald_Reci_Weight), intent(inout) :: this
    end subroutine Null_Ewald_Reci_Weight_destroy

    pure subroutine Null_Ewald_Reci_Weight_set(this)
        class(Null_Ewald_Reci_Weight), intent(inout) :: this
    end subroutine Null_Ewald_Reci_Weight_set

    pure function Null_Ewald_Reci_Weight_get_reci_numbers(this) result(reci_numbers)
        class(Null_Ewald_Reci_Weight), intent(in) :: this
        integer :: reci_numbers(num_dimensions)
        reci_numbers = 0
    end function Null_Ewald_Reci_Weight_get_reci_numbers

    pure real(DP) function Null_Ewald_Reci_Weight_get(this, n_1, n_2, n_3) result(weight)
        class(Null_Ewald_Reci_Weight), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3
        weight = 0._DP
    end function Null_Ewald_Reci_Weight_get

!end implementation Null_Ewald_Reci_Weight

end module class_ewald_reci_weight
