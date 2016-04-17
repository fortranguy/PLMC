module class_des_reci_weight

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use class_periodic_box, only: Abstract_Periodic_Box
use class_permittivity, only: Abstract_Permittivity
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use class_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter

implicit none

private

    type, abstract, public :: Abstract_DES_Reci_Weight
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        integer :: reci_numbers(num_dimensions) = 0
        real(DP) :: permittivity = 0._DP
        class(Abstract_DES_Convergence_Parameter), pointer :: alpha => null()
        real(DP), allocatable :: weight(:, :, :)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset => Abstract_set
        procedure :: get => Abstract_get
        procedure, private :: set => Abstract_set
    end type Abstract_DES_Reci_Weight

    type, extends(Abstract_DES_Reci_Weight), public :: Concrete_DES_Reci_Weight

    end type Concrete_DES_Reci_Weight

    type, extends(Abstract_DES_Reci_Weight), public :: Null_DES_Reci_Weight
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: reset => Null_set
        procedure :: get => Null_get
    end type Null_DES_Reci_Weight

contains

!implementation Abstract_DES_Reci_Weight

    subroutine Abstract_construct(this, periodic_box, reciprocal_lattice, permittivity, alpha)
        class(Abstract_DES_Reci_Weight), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_DES_Convergence_Parameter), target, intent(in) :: alpha

        this%periodic_box => periodic_box
        this%reci_numbers = reciprocal_lattice%get_numbers()
        this%permittivity = permittivity%get()
        this%alpha => alpha
        allocate(this%weight(0:this%reci_numbers(1), 0:this%reci_numbers(2), &
                             0:this%reci_numbers(3)))
        this%weight = 0._DP
        call this%set()
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_DES_Reci_Weight), intent(inout) :: this

        if (allocated(this%weight)) deallocate(this%weight)
        this%alpha => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    pure subroutine Abstract_set(this)
        class(Abstract_DES_Reci_Weight), intent(inout) :: this

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

            if (n_1**2 + n_2**2 + n_3**2 > this%reci_numbers(1)**2) cycle

            if (n_1 == 0 .and. n_2 == 0 .and. n_3 == 0) then
                this%weight(n_1, n_2, n_3) = 0._DP
            else
                wave_squared = dot_product(wave_vector, wave_vector)
                this%weight(n_1, n_2, n_3) = exp(-wave_squared/alpha**2/4._DP) / &
                    (2._DP*this%permittivity*product(box_size) * wave_squared)
            end if
        end do
        end do
        end do
    end subroutine Abstract_set

    !> \[
    !>      w_\alpha(\vec{k}) =
    !>          \begin{cases}
    !>              0   & \text{if } \vec{k} = \vec{0} \\
    !>              \frac{e^{-k^2/4\alpha^2}}{2\epsilon V k^2}  & \text{else}
    !>          \end{cases}
    !> \]
    pure real(DP) function Abstract_get(this, n_1, n_2, n_3) result(weight)
        class(Abstract_DES_Reci_Weight), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3

        weight = this%weight(abs(n_1), abs(n_2), abs(n_3))
    end function Abstract_get

!end implementation Abstract_DES_Reci_Weight

!implementation Null_DES_Reci_Weight

    subroutine Null_construct(this, periodic_box, reciprocal_lattice, permittivity, alpha)
        class(Null_DES_Reci_Weight), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_DES_Convergence_Parameter), target, intent(in) :: alpha
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DES_Reci_Weight), intent(inout) :: this
    end subroutine Null_destroy

    pure subroutine Null_set(this)
        class(Null_DES_Reci_Weight), intent(inout) :: this
    end subroutine Null_set

    pure real(DP) function Null_get(this, n_1, n_2, n_3) result(weight)
        class(Null_DES_Reci_Weight), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3
        weight = 0._DP
    end function Null_get

!end implementation Null_DES_Reci_Weight

end module class_des_reci_weight
