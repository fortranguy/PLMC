module classes_des_reci_weight

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use classes_box_size_memento, only: Abstract_Box_Size_Memento
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use classes_permittivity, only: Abstract_Permittivity
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter

implicit none

private

    type, abstract, public :: Abstract_DES_Reci_Weight
    private
        class(Abstract_Box_Size_Memento), pointer :: box_size_memento => null()
        integer :: reci_numbers(num_dimensions) = 0
        real(DP) :: permittivity = 0._DP
        real(DP) :: alpha_x_box_edge = 0._DP
        real(DP), allocatable :: weight(:, :, :)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: target => Abstract_target
        procedure :: reset => Abstract_reset
        procedure :: get => Abstract_get
    end type Abstract_DES_Reci_Weight

    type, extends(Abstract_DES_Reci_Weight), public :: Concrete_DES_Reci_Weight

    end type Concrete_DES_Reci_Weight

    type, extends(Abstract_DES_Reci_Weight), public :: Null_DES_Reci_Weight
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: target => Null_target
        procedure :: reset => Null_reset
        procedure :: get => Null_get
    end type Null_DES_Reci_Weight

contains

!implementation Abstract_DES_Reci_Weight

    !> @note Since box_size_memento may not be initialised,
    !> [[classes_des_reci_weight:Abstract_reset]] will be delayed.
    subroutine Abstract_construct(this, box_size_memento, reciprocal_lattice, permittivity, alpha)
        class(Abstract_DES_Reci_Weight), intent(out) :: this
        class(Abstract_Box_Size_Memento), intent(in) :: box_size_memento
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha

        call this%target(box_size_memento)
        this%reci_numbers = reciprocal_lattice%get_numbers()
        this%permittivity = permittivity%get()
        this%alpha_x_box_edge = alpha%get_times_box_edge()
        allocate(this%weight(0:this%reci_numbers(1), 0:this%reci_numbers(2), &
                             0:this%reci_numbers(3)))
        this%weight = 0._DP
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_DES_Reci_Weight), intent(inout) :: this

        if (allocated(this%weight)) deallocate(this%weight)
        this%box_size_memento => null()
    end subroutine Abstract_destroy

    subroutine Abstract_target(this, box_size_memento)
        class(Abstract_DES_Reci_Weight), intent(inout) :: this
        class(Abstract_Box_Size_Memento), target, intent(in) :: box_size_memento

        this%box_size_memento => box_size_memento
    end subroutine Abstract_target

    pure subroutine Abstract_reset(this)
        class(Abstract_DES_Reci_Weight), intent(inout) :: this

        integer :: n_1, n_2, n_3
        real(DP), dimension(num_dimensions) :: box_size, wave_vector
        real(DP) :: alpha, wave_squared

        box_size = this%box_size_memento%get()
        alpha = this%alpha_x_box_edge / box_size(1)
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
    end subroutine Abstract_reset

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

    subroutine Null_construct(this, box_size_memento, reciprocal_lattice, permittivity, alpha)
        class(Null_DES_Reci_Weight), intent(out) :: this
        class(Abstract_Box_Size_Memento), intent(in) :: box_size_memento
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DES_Reci_Weight), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_target(this, box_size_memento)
        class(Null_DES_Reci_Weight), intent(inout) :: this
        class(Abstract_Box_Size_Memento), target, intent(in) :: box_size_memento
    end subroutine Null_target

    pure subroutine Null_reset(this)
        class(Null_DES_Reci_Weight), intent(inout) :: this
    end subroutine Null_reset

    pure real(DP) function Null_get(this, n_1, n_2, n_3) result(weight)
        class(Null_DES_Reci_Weight), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3
        weight = 0._DP
    end function Null_get

!end implementation Null_DES_Reci_Weight

end module classes_des_reci_weight
